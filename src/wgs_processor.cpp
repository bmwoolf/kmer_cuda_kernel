#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include "fastq_parser.h"
#include <stdint.h>

extern "C" int run_kmer_kernel(const char* h_seq, uint32_t* h_out, int n, int k);

class WGSProcessor {
private:
    int k_;
    size_t batch_size_;
    std::vector<std::string> input_files_;
    
public:
    WGSProcessor(int k = 31, size_t batch_size = 1000000) 
        : k_(k), batch_size_(batch_size) {}
    
    void addInputFile(const std::string& filename) {
        input_files_.push_back(filename);
    }
    
    void processFiles() {
        std::cout << "=== WGS K-mer Processing ===" << std::endl;
        std::cout << "K-mer size: " << k_ << std::endl;
        std::cout << "Batch size: " << batch_size_ << " sequences" << std::endl;
        std::cout << "Input files: " << input_files_.size() << std::endl;
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        size_t total_sequences = 0;
        size_t total_kmers = 0;
        
        for (const auto& filename : input_files_) {
            std::cout << "\nProcessing: " << filename << std::endl;
            
            try {
                auto file_stats = processSingleFile(filename);
                total_sequences += file_stats.sequences_processed;
                total_kmers += file_stats.kmers_generated;
                
                std::cout << "  Sequences processed: " << file_stats.sequences_processed << std::endl;
                std::cout << "  K-mers generated: " << file_stats.kmers_generated << std::endl;
                std::cout << "  Processing time: " << std::fixed << std::setprecision(2) 
                          << file_stats.processing_time.count() << " seconds" << std::endl;
                
            } catch (const std::exception& e) {
                std::cerr << "Error processing " << filename << ": " << e.what() << std::endl;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\n=== Summary ===" << std::endl;
        std::cout << "Total sequences processed: " << total_sequences << std::endl;
        std::cout << "Total k-mers generated: " << total_kmers << std::endl;
        std::cout << "Total processing time: " << total_time.count() << " seconds" << std::endl;
        std::cout << "Average sequences per second: " 
                  << (total_time.count() > 0 ? total_sequences / total_time.count() : 0) << std::endl;
    }
    
private:
    struct FileStats {
        size_t sequences_processed;
        size_t kmers_generated;
        std::chrono::duration<double> processing_time;
    };
    
    FileStats processSingleFile(const std::string& filename) {
        FileStats stats = {0, 0, std::chrono::duration<double>(0)};
        auto start_time = std::chrono::high_resolution_clock::now();
        
        FastqParser parser(filename);
        std::string sequence, quality;
        std::vector<std::string> batch;
        batch.reserve(batch_size_);
        
        while (parser.readNextSequence(sequence, quality)) {
            std::string normalized = dna_utils::normalizeSequence(sequence);
            if (!normalized.empty() && normalized.length() >= k_) {
                batch.push_back(normalized);
                stats.sequences_processed++;
                
                if (batch.size() >= batch_size_) {
                    processBatch(batch, stats);
                    batch.clear();
                }
            }
        }
        
        // Process remaining sequences
        if (!batch.empty()) {
            processBatch(batch, stats);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        stats.processing_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
        
        return stats;
    }
    
    void processBatch(const std::vector<std::string>& sequences, FileStats& stats) {
        for (const auto& seq : sequences) {
            if (seq.length() >= k_) {
                // Process k-mers for this sequence
                std::vector<uint32_t> kmers(seq.length() - k_ + 1);
                
                int result = run_kmer_kernel(seq.c_str(), kmers.data(), seq.length(), k_);
                if (result == 0) {
                    stats.kmers_generated += kmers.size();
                    
                    // Optional: Print first few k-mers for verification
                    if (stats.sequences_processed <= 3) {
                        std::cout << "  Sample k-mers for sequence " << stats.sequences_processed << ":" << std::endl;
                        for (size_t i = 0; i < std::min(size_t(5), kmers.size()); ++i) {
                            std::cout << "    " << i << ": " << kmers[i] << std::endl;
                        }
                    }
                } else {
                    std::cerr << "CUDA kernel failed for sequence " << stats.sequences_processed << std::endl;
                }
            }
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <k-mer_size> [batch_size] [fastq_file1] [fastq_file2] ..." << std::endl;
        std::cout << "Example: " << argv[0] << " 31 1000000 ~/wgs_data/*.fastq.gz" << std::endl;
        return 1;
    }
    
    int k = std::stoi(argv[1]);
    size_t batch_size = (argc > 2) ? std::stoul(argv[2]) : 1000000;
    
    WGSProcessor processor(k, batch_size);
    
    // Add input files
    for (int i = 3; i < argc; ++i) {
        processor.addInputFile(argv[i]);
    }
    
    // If no files specified, use a test file
    if (argc <= 3) {
        std::cout << "No input files specified. Testing with a small sample..." << std::endl;
        
        // Create a test FASTQ file for demonstration
        std::string test_file = "test_sample.fastq";
        std::ofstream test_out(test_file);
        test_out << "@test_read_1\n";
        test_out << "ACGTACGTACGTACGTACGTACGTACGTACGT\n";
        test_out << "+\n";
        test_out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        test_out << "@test_read_2\n";
        test_out << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n";
        test_out << "+\n";
        test_out << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
        test_out.close();
        
        processor.addInputFile(test_file);
    }
    
    processor.processFiles();
    
    return 0;
} 