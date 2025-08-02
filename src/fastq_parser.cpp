#include "fastq_parser.h"
#include <algorithm>
#include <cctype>
#include <cmath>

FastqParser::FastqParser(const std::string& filename) : filename_(filename), is_gzipped_(false) {
    // Check if file is gzipped by looking at magic bytes
    std::ifstream test_file(filename, std::ios::binary);
    if (test_file.is_open()) {
        unsigned char magic[2];
        test_file.read(reinterpret_cast<char*>(magic), 2);
        is_gzipped_ = (magic[0] == 0x1f && magic[1] == 0x8b);
        test_file.close();
    }
    
    if (is_gzipped_) {
        gzfile_ = std::make_unique<gzFile>(gzopen(filename.c_str(), "r"));
        if (!*gzfile_) {
            throw std::runtime_error("Cannot open gzipped file: " + filename);
        }
    } else {
        file_ = std::make_unique<std::ifstream>(filename);
        if (!file_->is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
    }
}

FastqParser::~FastqParser() {
    if (is_gzipped_ && gzfile_) {
        gzclose(*gzfile_);
    }
}

bool FastqParser::readNextSequence(std::string& sequence, std::string& quality) {
    std::string line;
    bool success = false;
    
    // Read header line (starts with @)
    if (is_gzipped_) {
        char buffer[1024];
        if (gzgets(*gzfile_, buffer, sizeof(buffer)) == nullptr) {
            return false;
        }
        line = buffer;
    } else {
        if (!std::getline(*file_, line)) {
            return false;
        }
    }
    
    if (line.empty() || line[0] != '@') {
        return false;
    }
    
    // Read sequence line
    if (is_gzipped_) {
        char buffer[1024];
        if (gzgets(*gzfile_, buffer, sizeof(buffer)) == nullptr) {
            return false;
        }
        sequence = buffer;
    } else {
        if (!std::getline(*file_, sequence)) {
            return false;
        }
    }
    
    // Remove newline if present
    if (!sequence.empty() && sequence.back() == '\n') {
        sequence.pop_back();
    }
    
    // Read + line (optional)
    if (is_gzipped_) {
        char buffer[1024];
        if (gzgets(*gzfile_, buffer, sizeof(buffer)) == nullptr) {
            return false;
        }
        line = buffer;
    } else {
        if (!std::getline(*file_, line)) {
            return false;
        }
    }
    
    // Read quality line
    if (is_gzipped_) {
        char buffer[1024];
        if (gzgets(*gzfile_, buffer, sizeof(buffer)) == nullptr) {
            return false;
        }
        quality = buffer;
    } else {
        if (!std::getline(*file_, quality)) {
            return false;
        }
    }
    
    // Remove newline if present
    if (!quality.empty() && quality.back() == '\n') {
        quality.pop_back();
    }
    
    return true;
}

size_t FastqParser::getApproximateSequenceCount() {
    // This is a rough estimate based on file size
    // In practice, you'd want to do a full scan for accurate count
    std::ifstream file(filename_, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        return 0;
    }
    
    size_t file_size = file.tellg();
    file.close();
    
    // Rough estimate: assume average read length of 150bp + overhead
    // This is very approximate and depends on the actual data
    return file_size / 200; // Rough estimate
}

namespace dna_utils {

std::string normalizeSequence(const std::string& sequence) {
    std::string normalized;
    normalized.reserve(sequence.length());
    
    for (char c : sequence) {
        char upper = std::toupper(c);
        if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' || upper == 'N') {
            normalized += upper;
        }
        // Skip invalid characters
    }
    
    return normalized;
}

bool isValidDNA(const std::string& sequence) {
    for (char c : sequence) {
        char upper = std::toupper(c);
        if (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T' && upper != 'N') {
            return false;
        }
    }
    return true;
}

SequenceStats analyzeSequences(const std::string& filename) {
    SequenceStats stats = {0, 0, SIZE_MAX, 0, 0.0};
    
    try {
        FastqParser parser(filename);
        std::string sequence, quality;
        
        while (parser.readNextSequence(sequence, quality)) {
            std::string normalized = normalizeSequence(sequence);
            if (!normalized.empty()) {
                stats.total_sequences++;
                stats.total_bases += normalized.length();
                stats.min_length = std::min(stats.min_length, normalized.length());
                stats.max_length = std::max(stats.max_length, normalized.length());
            }
        }
        
        if (stats.total_sequences > 0) {
            stats.avg_length = static_cast<double>(stats.total_bases) / stats.total_sequences;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error analyzing sequences: " << e.what() << std::endl;
    }
    
    return stats;
}

} // namespace dna_utils 