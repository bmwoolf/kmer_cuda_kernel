// main CUDA kernel
#include <stdint.h>

__device__ __forceinline__ uint32_t encode_kmer(const char* seq, int k, int start) {
    uint32_t kmer = 0;
    for (int i = 0; i < k; ++i) {
        char base = seq[start + i];
        uint32_t val = 0;
        if (base == 'C') val = 1;
        else if (base == 'G') val = 2;
        else if (base == 'T') val = 3;
        kmer = (kmer << 2) | val;
    }
    return kmer;
}

__global__ void kmer_kernel(const char* seq, uint32_t* out, int n, int k) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i + k <= n) {
        out[i] = encode_kmer(seq, k, i);
    }
}

extern "C" void run_kmer_kernel(const char* h_seq, uint32_t* h_out, int n, int k) {
    char* d_seq;
    uint32_t* d_out;

    cudaMalloc(&d_seq, n * sizeof(char));
    cudaMalloc(&d_out, (n - k + 1) * sizeof(char), cudaMemcpyHostToDevice);

    cudaMemcpy(d_seq, h_seq, n * sizeof(char), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;
    kmer_kernel<<<numBlocks, blockSize>>>(d_seq, d_out, n, k);

    cudaMemcpy(h_out, d_out, (n - k + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost);

    cudaFree(d_seq);
    cudaFree(d_out);
}