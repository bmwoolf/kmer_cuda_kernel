// host code
#include <iostream>
#include <vector>
#include <string>
#include <stdint.h>

extern "C" int run_kmer_kernel(const char* h_seq, uint32_t* h_out, int n, int k);

int main() {
    std::string seq = "ACGTATCGGATCA";
    int k = 3;
    int n = seq.size();
    std::vector<uint32_t> out(n - k + 1);

    int result = run_kmer_kernel(seq.c_str(), out.data(), n, k);
    
    if (result != 0) {
        std::cerr << "CUDA kernel execution failed!" << std::endl;
        return 1;
    }

    for (int i = 0; i < n - k + 1; ++i) {
        std::cout << "K-mer[" << i << "]: " << out[i] << std::endl;
    }

    return 0;
}