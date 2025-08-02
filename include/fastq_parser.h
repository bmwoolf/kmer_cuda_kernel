#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <zlib.h>
#include <memory>

class FastqParser {
private:
    std::string filename_;
    std::unique_ptr<gzFile> gzfile_;
    std::unique_ptr<std::ifstream> file_;
    bool is_gzipped_;
    
public:
    FastqParser(const std::string& filename);
    ~FastqParser();
    
    // Read next sequence from FASTQ file
    bool readNextSequence(std::string& sequence, std::string& quality);
    
    // Get total number of sequences (approximate)
    size_t getApproximateSequenceCount();
    
    // Check if file is gzipped
    bool isGzipped() const { return is_gzipped_; }
    
    // Get filename
    const std::string& getFilename() const { return filename_; }
};

// Utility functions for DNA processing
namespace dna_utils {
    // Convert DNA sequence to uppercase and filter valid bases
    std::string normalizeSequence(const std::string& sequence);
    
    // Check if sequence contains valid DNA bases
    bool isValidDNA(const std::string& sequence);
    
    // Get sequence length statistics
    struct SequenceStats {
        size_t total_sequences;
        size_t total_bases;
        size_t min_length;
        size_t max_length;
        double avg_length;
    };
    
    SequenceStats analyzeSequences(const std::string& filename);
} 