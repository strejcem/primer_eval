# Primer Eval

A hyper-fast, memory-efficient 16S primer coverage evaluator written in Rust. 

Unlike traditional Python scripts that load massive databases into RAM or rely on slow regex engines, `primer_eval` uses bitwise IUPAC matching and streams data directly from disk. This allows it to evaluate gigabyte-sized FASTA/FASTQ databases against degenerate primers in seconds while keeping peak RAM usage nearly flat (< 20 MB).

> **Note:** This tool was *vibecoded* (built with AI assistance), but the underlying logic, edge-cases, and performance metrics have been thoroughly tested and validated.

## Installation

### Option 1: Download the Pre-compiled Binary
The easiest way to use `primer_eval` is to grab the pre-compiled executable from the **Releases** tab on the right side of this GitHub page. Just download it, make it executable (`chmod +x primer_eval`), and run it!

### Option 2: Build from Source
If you prefer to compile it yourself, you will need [Rust and Cargo](https://rustup.rs/) installed.

1. Clone the repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/primer_eval.git](https://github.com/YOUR_USERNAME/primer_eval.git)
   cd primer_eval
   ```
2. Build the optimized release binary:
   ```bash
   cargo build --release
   ```
3. The executable will be located at `target/release/primer_eval`.

## Usage

Run the tool via the command line. By default, it will rapidly evaluate the coverage percentage and output a summary. 

```bash
./primer_eval -p GTGCCAGCMGCCGCGGTAA -d silva_db.fasta -m 2 -s 5
```

### Generating Detailed Hit Locations
If you want to know exactly *where* the primer hit in every sequence, use the `-o` or `--output` flag. This will generate a TSV file containing the Sequence ID, Strand (+/-), 1-based Position, and the number of mismatches (MM).

```bash
./primer_eval -p GTGCCAGCMGCCGCGGTAA -d silva_db.fasta -m 2 -s 5 -o hits.tsv
```

### Arguments

| Flag | Name | Default | Description |
|---|---|---|---|
| `-p` | `--primer` | **Required** | The primer sequence (e.g., `GTGCCAGCMGCCGCGGTAA`). Accepts IUPAC degenerate bases. |
| `-d` | `--database` | **Required** | Path to the target FASTA/FASTQ database file. |
| `-o` | `--output` | None | Optional path to save a detailed hit locations TSV file. |
| `-m` | `--mismatches`| `2` | Max mismatches allowed in the 5' fuzzy region. |
| `-s` | `--strict3` | `5` | Length of the strict 0-mismatch 3' anchor. |
| `-t` | `--threads` | All Cores | Number of CPU cores to allocate. |

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
