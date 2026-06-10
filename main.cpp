#include <iostream>
#include "fasta.h"
#include <string>
#include <fstream>
using namespace std;

//Struct to hold each record
void printStats(const ProteinStats& stats) {
	cout << "Length:            " << stats.length << endl;
	cout << "Molecular  Weight: " << stats.molecularWeight << " Da" << endl;
	cout << "Composition:" endl;
	for (const auto& pair : stats.aminoAcidComposition) {
		cout << "  " << pair.first << ": " << pair.second << endl;
	}
}

int main(){
/*
 * This analyzer has five goals.
 * 1- Import and Parse a FASTA File (and the codon table), and extract sequences
 * 2- Compute basic stats: Length, Amino Acid composition, molecular weight
 * 3- Search for motifs using a trie
 * 4- Translate a DNA sequence to a protein sequence using the codon table
 * 5- Write Smith-Waterman local alignment (This will have to come later, of course)
 

*/
	//Parse File
	vector<FastaRecord> records = parseFasta("EXAMPLE_1.fasta");
	

	// Print parsed records:
	for (const auto& rec : records) {
		cout << "Header:  " << rec.header << endl;
		cout << "Sequence " << rec.sequence << endl;
		cout << "Length   " << rec.sequence.size() << "bp" << endl;
	}
	return 0;
	
}
