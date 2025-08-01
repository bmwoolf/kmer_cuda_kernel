// main CUDA kernel
#include <stdint.h>
#include <cuda_runtime.h>
#include <iostream>

__device__ __forceinline__ uint32_t encode_kmer(const char* seq, int k, int start) {
    uint32_t kmer = 0;
    for (int i = 0; i < k; ++i) {
        char base = seq[start + i];
        uint32_t val = 0;
        if (base == 'A') val = 0;
        else if (base == 'C') val = 1;
        else if (base == 'G') val = 2;
        else if (base == 'T') val = 3;
        // For any other character (N, etc.), treat as A
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

extern "C" int run_kmer_kernel(const char* h_seq, uint32_t* h_out, int n, int k) {
    char* d_seq;
    uint32_t* d_out;
    cudaError_t err;

    // Allocate device memory
    err = cudaMalloc(&d_seq, n * sizeof(char));
    if (err != cudaSuccess) {
        std::cerr << "CUDA malloc error for seq: " << cudaGetErrorString(err) << std::endl;
        return -1;
    }

    err = cudaMalloc(&d_out, (n - k + 1) * sizeof(uint32_t));
    if (err != cudaSuccess) {
        std::cerr << "CUDA malloc error for out: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_seq);
        return -1;
    }

    // Copy data to device
    err = cudaMemcpy(d_seq, h_seq, n * sizeof(char), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        std::cerr << "CUDA memcpy error: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_seq);
        cudaFree(d_out);
        return -1;
    }

    // Launch kernel
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;
    kmer_kernel<<<numBlocks, blockSize>>>(d_seq, d_out, n, k);

    // Check for kernel launch errors
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "Kernel launch error: " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_seq);
        cudaFree(d_out);
        return -1;
    }

    // Copy results back to host
    err = cudaMemcpy(h_out, d_out, (n - k + 1) * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        std::cerr << "CUDA memcpy error (device to host): " << cudaGetErrorString(err) << std::endl;
        cudaFree(d_seq);
        cudaFree(d_out);
        return -1;
    }

    // Cleanup
    cudaFree(d_seq);
    cudaFree(d_out);
    
    return 0;
}