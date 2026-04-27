use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::sync::mpsc;
use std::thread::{self, available_parallelism};
use std::time::Instant;

/// Calculate fast 16S primer coverage with strict 3' anchoring.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// The primer sequence to evaluate (e.g., GTGCCAGCMGCCGCGGTAA).
    /// Accepts standard IUPAC degenerate bases.
    #[arg(short, long)]
    primer: String,

    /// Path to the target FASTA/FASTQ database file to search against.
    #[arg(short, long)]
    database: String,

    /// Path to save the detailed hit locations TSV file (Optional). 
    /// Outputs every sequence. Unmatched sequences are padded with NA.
    /// 'Position' in the output is 1-based (the first base of the sequence is index 1).
    #[arg(short, long)]
    output: Option<String>,

    /// Maximum number of mismatches permitted in the 5' fuzzy region of the primer.
    #[arg(short, long, default_value_t = 2)]
    mismatches: usize,

    /// Length of the strict 0-mismatch 3' anchor region.
    /// Primer matches will fail instantly if there are any mismatches in this region.
    #[arg(short, long, default_value_t = 5)]
    strict3: usize,

    /// Number of CPU cores to use for multiprocessing. 
    /// Defaults to all available physical/logical cores on the machine.
    #[arg(short, long)]
    threads: Option<usize>,
}

// Convert IUPAC nucleotide characters to bitmasks
const fn iupac_to_bitmask(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'G' | b'g' => 4,
        b'T' | b't' | b'U' | b'u' => 8,
        b'R' | b'r' => 1 | 4,
        b'Y' | b'y' => 2 | 8,
        b'S' | b's' => 2 | 4,
        b'W' | b'w' => 1 | 8,
        b'K' | b'k' => 4 | 8,
        b'M' | b'm' => 1 | 2,
        b'B' | b'b' => 2 | 4 | 8,
        b'D' | b'd' => 1 | 4 | 8,
        b'H' | b'h' => 1 | 2 | 8,
        b'V' | b'v' => 1 | 2 | 4,
        b'N' | b'n' | b'I' | b'i' => 15,
        _ => 0,
    }
}

// Reverse complement a bitmask
fn rc_bitmask(b: u8) -> u8 {
    let mut rc = 0;
    if b & 1 != 0 { rc |= 8; } // A -> T
    if b & 2 != 0 { rc |= 4; } // C -> G
    if b & 4 != 0 { rc |= 2; } // G -> C
    if b & 8 != 0 { rc |= 1; } // T -> A
    rc
}

const fn build_lut() -> [u8; 256] {
    let mut lut = [0; 256];
    let mut i = 0;
    while i < 256 {
        lut[i] = iupac_to_bitmask(i as u8);
        i += 1;
    }
    lut
}
static SEQ_LUT: [u8; 256] = build_lut();

struct PrimerMatcher {
    fwd_masks: Vec<u8>,
    rev_masks: Vec<u8>,
    l5: usize,
    l3: usize,
    max_mismatches: usize,
}

impl PrimerMatcher {
    fn new(primer: &str, strict3: usize, max_mismatches: usize) -> Self {
        let p_bytes = primer.as_bytes();
        let l5 = p_bytes.len().saturating_sub(strict3);
        let l3 = strict3;

        let fwd_masks: Vec<u8> = p_bytes.iter().map(|&b| iupac_to_bitmask(b)).collect();
        let rev_masks: Vec<u8> = p_bytes.iter().rev().map(|&b| rc_bitmask(iupac_to_bitmask(b))).collect();

        Self { fwd_masks, rev_masks, l5, l3, max_mismatches }
    }

    /// Fast path: Returns true the instant it finds a match.
    fn is_hit(&self, seq: &[u8]) -> bool {
        let len = self.fwd_masks.len();
        if seq.len() < len { return false; }
        for pos in 0..=(seq.len() - len) {
            if self.check_fwd(seq, pos).is_some() || self.check_rev(seq, pos).is_some() {
                return true;
            }
        }
        false
    }

    /// Detailed path: Scans the entire sequence and returns hits: (Position, Strand, Mismatches)
    fn find_hits(&self, seq: &[u8]) -> Vec<(usize, char, usize)> {
        let mut hits = Vec::new();
        let len = self.fwd_masks.len();
        if seq.len() < len { return hits; }
        
        for pos in 0..=(seq.len() - len) {
            if let Some(mm) = self.check_fwd(seq, pos) {
                hits.push((pos + 1, '+', mm)); // 1-based start position
            } else if let Some(mm) = self.check_rev(seq, pos) {
                hits.push((pos + 1, '-', mm)); // 1-based start position
            }
        }
        hits
    }

    #[inline(always)]
    fn check_fwd(&self, seq: &[u8], pos: usize) -> Option<usize> {
        // Fast fail on strict 3' end
        for i in self.l5..self.fwd_masks.len() {
            if (SEQ_LUT[seq[pos + i] as usize] & self.fwd_masks[i]) == 0 { return None; }
        }
        let mut mismatches = 0;
        for i in 0..self.l5 {
            if (SEQ_LUT[seq[pos + i] as usize] & self.fwd_masks[i]) == 0 {
                mismatches += 1;
                if mismatches > self.max_mismatches { return None; }
            }
        }
        Some(mismatches)
    }

    #[inline(always)]
    fn check_rev(&self, seq: &[u8], pos: usize) -> Option<usize> {
        // Reverse complement flips primer sides. Strict anchor is now at the 5' end of the slice.
        for i in 0..self.l3 {
            if (SEQ_LUT[seq[pos + i] as usize] & self.rev_masks[i]) == 0 { return None; }
        }
        let mut mismatches = 0;
        for i in self.l3..self.rev_masks.len() {
            if (SEQ_LUT[seq[pos + i] as usize] & self.rev_masks[i]) == 0 {
                mismatches += 1;
                if mismatches > self.max_mismatches { return None; }
            }
        }
        Some(mismatches)
    }
}

fn main() {
    let args = Args::parse();

    if args.primer.len() <= args.strict3 {
        eprintln!("Error: Primer length ({}) must be strictly greater than the 3' anchor length ({}).", 
                  args.primer.len(), args.strict3);
        std::process::exit(1);
    }

    let max_cores = available_parallelism().map(|n| n.get()).unwrap_or(1);
    let threads = args.threads.unwrap_or(max_cores).clamp(1, max_cores);

    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

    println!("--- Initialization ---");
    println!("Primer: {}", args.primer);
    println!("Database: {}", args.database);
    println!("5' Mismatches allowed: {}", args.mismatches);
    println!("3' Strict anchor length: {}", args.strict3);
    println!("Threads allocated: {}", threads);
    if args.output.is_some() {
        println!("Detailed output: Enabled (Requires full sequence scans)");
    } else {
        println!("Detailed output: Disabled (Using fast short-circuit evaluation)");
    }
    
    let matcher = PrimerMatcher::new(&args.primer, args.strict3, args.mismatches);
    let write_locations = args.output.is_some();

    // Setup an I/O background thread for writing detailed locations TSV
    let mut detail_writer = None;
    if let Some(ref out_file) = args.output {
        let (tx, rx) = mpsc::channel::<String>();
        let out_file_clone = out_file.clone();
        let handle = thread::spawn(move || {
            let mut file = File::create(&out_file_clone).expect("Could not create output file");
            // R-friendly and Python-friendly header
            writeln!(file, "Sequence_ID\tStrand\tPosition\tMismatches").unwrap();
            for msg in rx {
                file.write_all(msg.as_bytes()).unwrap();
            }
        });
        detail_writer = Some((tx, handle));
    }

    println!("\nStreaming database and evaluating sequences...");
    let start_time = Instant::now();

    let mut reader = match parse_fastx_file(&args.database) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error: Could not open database file '{}': {}", args.database, e);
            std::process::exit(1);
        }
    };

    let mut total_seqs: usize = 0;
    let mut total_hits: usize = 0;
    let chunk_size = 100_000; 
    let mut sequence_chunk = Vec::with_capacity(chunk_size);

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid sequence record");
        // Store both the ID and the sequence
        sequence_chunk.push((seqrec.id().to_vec(), seqrec.normalize(false).into_owned()));

        if sequence_chunk.len() >= chunk_size {
            let (chunk_hits, out_text) = process_chunk(&sequence_chunk, &matcher, write_locations);
            total_hits += chunk_hits;
            total_seqs += sequence_chunk.len();
            
            if let (Some(text), Some((tx, _))) = (out_text, detail_writer.as_ref()) {
                let _ = tx.send(text);
            }
            sequence_chunk.clear();
        }
    }

    // Process remainder
    if !sequence_chunk.is_empty() {
        let (chunk_hits, out_text) = process_chunk(&sequence_chunk, &matcher, write_locations);
        total_hits += chunk_hits;
        total_seqs += sequence_chunk.len();
        
        if let (Some(text), Some((tx, _))) = (out_text, detail_writer.as_ref()) {
            let _ = tx.send(text);
        }
    }

    // Safely shut down the background writer thread
    if let Some((tx, handle)) = detail_writer {
        drop(tx); // Close the channel
        handle.join().expect("Writer thread panicked");
        println!("Detailed output successfully saved to {}", args.output.unwrap());
    }

    if total_seqs == 0 {
        eprintln!("Error: No sequences found in the provided database.");
        std::process::exit(1);
    }

    let elapsed = start_time.elapsed().as_secs_f64();
    let coverage = (total_hits as f64 / total_seqs as f64) * 100.0;

    let result_text = format!(
        "--- Coverage Results ---\n\
        Primer: {}\n\
        Database: {}\n\
        Total Sequences: {}\n\
        Sequences with Hit: {}\n\
        Coverage: {:.2}%\n\
        Processed in: {:.2} seconds\n",
        args.primer, args.database, total_seqs, total_hits, coverage, elapsed
    );

    println!("\n{}", result_text);
}

fn process_chunk(
    chunk: &[(Vec<u8>, Vec<u8>)], 
    matcher: &PrimerMatcher, 
    write_locations: bool
) -> (usize, Option<String>) {
    if write_locations {
        // Detailed path: Extract hits or pad with NA
        let results: Vec<(usize, String)> = chunk.par_iter().map(|(id_bytes, seq)| {
            let hits = matcher.find_hits(seq);
            let id_str = String::from_utf8_lossy(id_bytes);
            
            if hits.is_empty() { 
                // Return 0 hits, and an NA padded row
                return (0, format!("{}\tNA\tNA\tNA\n", id_str)); 
            }
            
            let mut lines = String::new();
            for (pos, strand, mm) in hits {
                lines.push_str(&format!("{}\t{}\t{}\t{}\n", id_str, strand, pos, mm));
            }
            // If there's at least one hit, this sequence counts towards coverage
            (1, lines)
        }).collect();
        
        // Sum the hit counts
        let hit_count: usize = results.iter().map(|(count, _)| count).sum();
        
        // Combine the strings
        let mut combined_text = String::new();
        for (_, text) in results { 
            combined_text.push_str(&text); 
        }
        
        (hit_count, Some(combined_text))
    } else {
        // Fast path: Just short-circuit count hits
        let hit_count: usize = chunk.par_iter().map(|(_, seq)| {
            if matcher.is_hit(seq) { 1 } else { 0 }
        }).sum();
        (hit_count, None)
    }
}