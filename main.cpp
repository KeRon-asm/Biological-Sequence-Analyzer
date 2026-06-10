#include <iostream>
#include "fasta.h"
#include <string>
#include <fstream>
#include <mutex>
#include <thread>
#include <functional>
#include <vector>
using namespace std;

//Struct to hold each record
void printStats(const ProteinStats& stats) {
	cout << "Length:            " << stats.length << endl;
	cout << "Molecular  Weight: " << stats.molecularWeight << " Da" << endl;
	cout << "Composition:" << endl;
	for (const auto& pair : stats.aminoAcidComposition) {
		// v
		cout << "  " << pair.first << ": " << pair.second << endl;
	}
}

int main(){
'''
This will be commented out once Im sure how I want to implement the tool:
Allows adding an argument instead of using the example file provided

int main(int  argc, char* argv[]) {
	if (argc < 2) {
		cerr << "Error, try the following output: ./analyzer <filename.fasta>" << endl;
		return 1;
	}
}

	string filename = argv[1];
	vector<FastaRecord> records = parseFasta(filename);
	
'''



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
	
	vector<thread> threads;
	for (const auto& rec : records) {
		threads.emplace_back(analyzeRecord, ref(rec));
	}
	for (auto& t : threads) {
		t.join();
	}

	// Print parsed records:
	for (const auto& rec : records) {
		cout << "Header:  " << rec.header << endl;
		cout << "Sequence " << rec.sequence << endl;
		cout << "Length   " << rec.sequence.size() << "bp" << endl;
	}
	return 0;
	
}
