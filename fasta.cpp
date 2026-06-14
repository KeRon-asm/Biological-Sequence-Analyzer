#include "fasta.h"
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
using namespace std;
vector<FastaRecord> parseFasta(const string& filename) {
    ifstream iFile(filename);
    if (!iFile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return {};
    }
    // Begin parsing file
	vector<FastaRecord> records;
	string line;
	FastaRecord current;
	//!!!
	while (getline(iFile, line)){
		if (line.empty()) continue; // Just skip any lines that are blank

		if (line[0] == '>') { // '>' character signals the header
			// save previous record if there is one
			if (!current.header.empty()){
				records.push_back(current);
				current = FastaRecord{}; //reset
			}
			current.header = line.substr(1);
		} else {
            current.sequence += line;
        }
	}
    if (!current.header.empty())
        records.push_back(current);
    return records;
}
string transcribeDNA(const string& dna) {
	// takes a string and swaps out the T for a U
	string rna = dna;
	for (char& c : rna) {
		// strings are mutable in C++
		if (c == 'T') c = 'U';
	}
	return rna;
}

string translateRNA(const string& rna){
	//take the rna string
	string composition;
	
	for (int i = 0; i+3 <= rna.size(); i+=3){
	string current_acid = rna.substr(i,3); 
	string codon;
	if (CODONTABLE.count(current_acid)) {
		codon =  CODONTABLE.at(current_acid);
	} else {
		codon = "UNKNOWN";
		cout << "Unknown Codon detected." << endl;
	}
	if (codon == "STOP") break;
	composition += codon;
	}
	return composition;
}
// Molecular weights of each amino acid in Daltons
inline const std::unordered_map<std::string, double> AA_WEIGHTS = {
    {"A", 89.09}, {"R", 174.20}, {"N", 132.12}, {"D", 133.10},
    {"C", 121.16}, {"Q", 146.15}, {"E", 147.13}, {"G", 75.03},
    {"H", 155.16}, {"I", 131.17}, {"L", 131.17}, {"K", 146.19},
    {"M", 149.21}, {"F", 165.19}, {"P", 115.13}, {"S", 105.09},
    {"T", 119.12}, {"W", 204.23}, {"Y", 181.19}, {"V", 117.15}
};

ProteinStats computeStats(const string& protein) {
    ProteinStats stats;
    
    stats.length = protein.size();

    // Goal: walk the protein string and count each amino acid
    //       into stats.aminoAcidComposition

    for (char c : protein) {
        string aminoAcid(1, c); // convert char to a 1-character string

        stats.aminoAcidComposition[aminoAcid]++;

        if (AA_WEIGHTS.count(aminoAcid)) {
            stats.molecularWeight += AA_WEIGHTS.at(aminoAcid);
        }
    }
        // Don't need to account for STOP codon, translateRNA already breaks out of the loop
    
    return stats;
    // Goal: use AA_WEIGHTS to accumulate stats.molecularWeight
}
std::mutex printMutex;
void analyzeRecord(const FastaRecord& rec) {
    // full pipeline
    // lock mutex before printing
}

//Smith Waterman Scoring constants

const int MATCH_SCORE = 2;
const int MISMATCH_SCORE = -1;
const int GAP_PENALTY = -1;

AlignmentResult smithWaterman(const string& seqA, const string& seqB) {
    int rows = seqA.size() + 1;
    int cols = seqB.size() + 1;

    AlignmentResult result;

    return result;
}

TrieNode* buildTrie(const vector<string>& motifs){
	TrieNode* root = new TrieNode();
	
	// for each motif in motifs:
	// start at root
	// for each character in motif:
	// ...

	return root
}
// Take a list of motifs, insert each one into the trie character by character
vector<string> searchMotifs(const string& protein, TrieNode* root);
// traverse through protein string, and at every starting position, follow ther trie as far as possible, recording any motifs

