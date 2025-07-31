# K-mer CUDA Kernel

This project aims to build a k-mer CUDA kernel from scratch. The core algorithm:
1. looks at a window of `k` letters in the input DNA sequence 
2. turns each letter into 2 bits
- A: 00
- C: 01
- T: 10
- G: 11
3. slides that window one base forward and repeats
4. stores that result as a number

Once this is done, you can use Smith-Waterman and the other algorithms on top of the converted sequence.

Our specific implementation will take in a 55GB series of FASTQ files from a whole-genome sequencer. 