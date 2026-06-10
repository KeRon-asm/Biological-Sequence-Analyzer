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
inline const unordered_map<string, double> AA_WEIGHTS = {
    {"ALA (A)", 89.09},
    {"ARG (R)", 174.20},
    {"ASN (N)", 132.12},
    {"ASP (D)", 133.10},
    {"CYS (C)", 121.16},
    {"GLN (Q)", 146.15},
    {"GLU (E)", 147.13},
    {"GLY (G)", 75.03},
    {"HIS (H)", 155.16},
    {"ILE (I)", 131.17},
    {"LEU (L)", 131.17},
    {"LYS (K)", 146.19},
    {"MET (M)", 149.21},
    {"PHE (F)", 165.19},
    {"PRO (P)", 115.13},
    {"SER (S)", 105.09},
    {"THR (T)", 119.12},
    {"TRP (W)", 204.23},
    {"TYR (Y)", 181.19},
    {"VAL (V)", 117.15}
};

ProteinStats computeStats(const string& protein) {
    ProteinStats stats;
    
    stats.length = protein.size();

    // Goal: walk the protein string and count each amino acid
    //       into stats.aminoAcidComposition

    for (int i = 0; i+7 <= protein.size(); i+= 7){
        string aminoAcid = protein.substr(i,7);

        // count composition
        stats.aminoAcidComposition[aminoAcid]++;

        // Accumulate molecular weight
        if (AA_WEIGHTS.count(aminoAcid)) {
            stats.molecularWeight += AA_WEIGHTS.at(aminoAcid);
        }
        // Don't need to account for STOP codon, translateRNA already breaks out of the loop
    }

    return stats;
    // Goal: use AA_WEIGHTS to accumulate stats.molecularWeight
    
    return stats;
}
std::mutex printMutex;
void analyzeRecord(const FastaRecord& rec) {
    // full pipeline
    // lock mutex before printing
}
