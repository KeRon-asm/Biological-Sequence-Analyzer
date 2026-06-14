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
inline const std::unordered_map<std::string, std::string> CODONTABLE = {
    {"AUU", "I"}, {"AUC", "I"}, {"AUA", "I"}, {"AUG", "M"},
    {"ACU", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},
    {"AAU", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"},
    {"AGU", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
    {"GUU", "V"}, {"GUC", "V"}, {"GUA", "V"}, {"GUG", "V"},
    {"GCU", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},
    {"GAU", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"},
    {"GGU", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"},
    {"UUU", "F"}, {"UUC", "F"}, {"UUA", "L"}, {"UUG", "L"},
    {"UCU", "S"}, {"UCC", "S"}, {"UCA", "S"}, {"UCG", "S"},
    {"UAU", "Y"}, {"UAC", "Y"}, {"UAA", "STOP"}, {"UAG", "STOP"},
    {"UGU", "C"}, {"UGC", "C"}, {"UGA", "STOP"}, {"UGG", "W"},
    {"CUU", "L"}, {"CUC", "L"}, {"CUA", "L"}, {"CUG", "L"},
    {"CCU", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},
    {"CAU", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
    {"CGU", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},
};

// Calculate Protein Stats
struct ProteinStats {
    int length;
    std::unordered_map<std::string, int> aminoAcidComposition;
    double molecularWeight;
};

ProteinStats computeStats(const std::string& protein);

	// Convert DNA to RNA
std::string transcribeDNA(const std::string& sequence);

std::string translateRNA(const std::string& sequence);

void analyzeRecord(const FastaRecord& rec);

struct AlignmentResult {
	std::string alignedA;
	std::string alignedB;
	int score;
};

AlignmentResult smithWaterman(const std::string& seqA, const std::string& seqB);

struct TrieNode {
    std::unordered_map<char, TrieNode*> children;
    bool isEndOfWord = false;
};
