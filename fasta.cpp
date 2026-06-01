#include "fasta.h"
#include <fstream>
#include <iostream>
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