#pragma once
#include <string>
#include <vector>
#include <unordered_map>

// should not use namespace here because it forces that choice on any file including this one
 
struct FastaRecord {
    std::string header;
    std::string sequence;
};

// Parse the file and return the record vectors
std::vector<FastaRecord> parseFasta(const std::string& filename);

// Build and return codon table
std::unordered_map<std::string, std::string> buildCodonTable();