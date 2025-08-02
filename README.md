# CUDA K-mer Kernel

A high-performance CUDA implementation for k-mer generation from DNA sequences, optimized for whole-genome sequencing (WGS) data processing.


## Core Algorithm

The k-mer kernel implements the following algorithm:

1. DNA Encoding: convert each base to 2 bits
   - A: 00, C: 01, G: 10, T: 11
2. Sliding Window: process k-length windows across the sequence
3. Parallel Processing: each GPU thread handles a different k-mer position
4. Bit Packing: store k-mers as 32-bit integers

## Building

### Prerequisites
- CUDA Toolkit (tested with CUDA 12.9)
- NVIDIA GPU with compute capability 5.0+
- GCC/G++ compiler
- zlib development library

### Installation
```bash
# Install dependencies (Ubuntu/Debian)
sudo apt install nvidia-cuda-toolkit g++ zlib1g-dev

# Build the project
make all
```

### Basic K-mer Generation
```bash
# Generate k-mers from a simple sequence
./kmer
```

### WGS Data Processing
```bash
# Process WGS FASTQ files
./wgs_processor 31 1000000 /path/to/wgs_data/*.fastq.gz
```

### Parameters
- `k-mer_size`: Length of k-mers (typically 31 for WGS)
- `batch_size`: Number of sequences to process per batch
- `fastq_files`: Input FASTQ files (supports .gz compression)


## Example Output

```
=== WGS K-mer Processing ===
K-mer size: 31
Batch size: 1000000 sequences
Input files: 1

Processing: /path/to/wgs_data.fastq.gz
  ✓ Sequences processed: 1000000
  ✓ K-mers generated: 50000000
  ✓ Processing time: 45.23 seconds
  ✓ Sequences/sec: 22109.4

=== Summary ===
✓ Total sequences processed: 1000000
✓ Total k-mers generated: 50000000
✓ Total processing time: 45 seconds
✓ Average sequences per second: 22109
```

### Testing
```bash
# Build and test basic functionality
make kmer && ./kmer

# Test with WGS data
make wgs_processor && ./wgs_processor 31 1000 test.fastq.gz
```

