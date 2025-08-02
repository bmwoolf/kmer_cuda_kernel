# CUDA K-mer Kernel Makefile
CC = nvcc
CFLAGS = -lcudart -lz
INCLUDES = -Iinclude

# Core k-mer kernel
kmer: examples/main.cpp src/kmer.cu
	$(CC) -o kmer examples/main.cpp src/kmer.cu $(CFLAGS)

# WGS processor (optional)
wgs_processor: src/wgs_processor.cpp src/fastq_parser.cpp src/kmer.cu
	$(CC) -o wgs_processor src/wgs_processor.cpp src/fastq_parser.cpp src/kmer.cu $(CFLAGS) $(INCLUDES)

# Build everything
all: kmer wgs_processor

# Clean build artifacts
clean:
	rm -f kmer wgs_processor

.PHONY: all clean