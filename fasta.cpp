#include "fasta.h"
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
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
	
	for (size_t i = 0; i+3 <= rna.size(); i+=3){
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
string detectSequenceType(const string& sequence) {
    // RNA is unambiguous — U only appears in RNA
    if (sequence.find('U') != string::npos) return "RNA";

    // Count characters that only appear in protein sequences
    // (B, D, E, F, H, I, J, K, L, M, N, O, P, Q, R, S, V, W, X, Y, Z)
    int proteinOnlyChars = 0;
    for (char c : sequence) {
        if (c=='E'||c=='F'||c=='I'||c=='L'||c=='M'||
            c=='P'||c=='Q'||c=='R'||c=='W'||c=='Y'||c=='D'||c=='H'||c=='K'||c=='N'||c=='S'||c=='V') {
            proteinOnlyChars++;
        }
    }

    // If more than 5% of characters are protein-only, it's a protein
    double ratio = (double)proteinOnlyChars / sequence.size();
    if (ratio > 0.05) return "PROTEIN";

    return "DNA";
}

void analyzeRecord(const FastaRecord& rec, 
                   const vector<string>& motifs, 
                   const string& querySequence) {    
    string type = detectSequenceType(rec.sequence);
    string protein;

    if (type == "DNA") {
        string rna = transcribeDNA(rec.sequence);
        protein = translateRNA(rna);
    } else if (type == "RNA") {
        protein = translateRNA(rec.sequence);
    } else {
        protein = rec.sequence;
    }

    ProteinStats stats = computeStats(protein);

    // Motif search
    
    TrieNode* root = buildTrie(motifs);
    vector<string> found = searchMotifs(protein, root);
    deleteTrie(root);

    lock_guard<mutex> lock(printMutex);
    cout << "\nHeader:           " << rec.header << endl;
    cout << "Sequence Type:    " << type << endl;
    cout << "Length:           " << stats.length << endl;
    cout << "Molecular Weight: " << stats.molecularWeight << " Da" << endl;
    cout << "Composition:" << endl;
    for (const auto& pair : stats.aminoAcidComposition) {
        cout << "  " << pair.first << ": " << pair.second << endl;
    }
    cout << "Motifs Found:" << endl;
    if (found.empty()) {
        cout << "  None" << endl;
    } else {
        for (const string& m : found) {
            cout << "  " << m << endl;
        }
    }
    if (!querySequence.empty()) {
        AlignmentResult alignment = smithWaterman(protein, querySequence);
        cout << "Alignment Score:  " << alignment.score << endl;
        cout << "Aligned Sequence: " << alignment.alignedA << endl;
        cout << "Query Aligned:    " << alignment.alignedB << endl;
    }
}

//Smith Waterman Scoring constants

const int MATCH_SCORE = 2;
const int MISMATCH_SCORE = -1;
const int GAP_PENALTY = -1;

AlignmentResult smithWaterman(const string& seqA, const string& seqB) {
    int rows = seqA.size() + 1;
    int cols = seqB.size() + 1;
    //Create matrix, initialized to 0
    vector<vector<int>> matrix(rows, vector<int>(cols,0));

    int bestScore = 0;
    int bestRow = 0;
    int bestCol = 0;

    //Fill matrix
    for (int i = 1; i < rows; i++){
        for (int j = 1; j < cols; j++) {

            // Calculate diagonal score
            int diagonal_score = 0;
            int up_score = 0;
            int left_score = 0;
            if (seqA[i-1] == seqB[j-1]) {
                diagonal_score=matrix[i-1][j-1] + MATCH_SCORE;
            } else {
                diagonal_score=matrix[i-1][j-1] + MISMATCH_SCORE;
            }
            // Up score
            up_score=matrix[i-1][j] + GAP_PENALTY;

            // Left score
            left_score=matrix[i][j-1] + GAP_PENALTY;

            matrix[i][j] = max(max(0,diagonal_score),max(up_score,left_score));
            if (bestScore<matrix[i][j]) {
                bestScore = matrix[i][j];
                bestRow = i;
                bestCol = j;
            }
        }
    }
    //Traceback
    string alignedA = "";
    string alignedB = "";
    int i = bestRow;
    int j = bestCol;
    while (i > 0 && j > 0 && matrix[i][j] != 0) {
        if ((seqA[i-1] == seqB[j-1]) && (matrix[i][j] == matrix [i-1][j-1] + MATCH_SCORE)){
            //If characters match and the score is consistent with the diagonal direction -> ensures we don't get accidental matches from up or left instead of diagonal
            alignedA.insert(0,1,seqA[i-1]); 
            alignedB.insert(0,1,seqB[j-1]);
            i--;
            j--;
        } else if (matrix[i][j] == matrix[i-1][j] + GAP_PENALTY) {
            //Match came from up_score, gap in seqB
            alignedA.insert(0,1,seqA[i-1]); 
            alignedB.insert(0,1,'-');
            i--;
        } else {
            //Match came from left_score, gap in seqA
            alignedA.insert(0,1,'-'); 
            alignedB.insert(0,1,seqB[j-1]);
            j--;
        }
    }

    AlignmentResult result;
    result.alignedA = alignedA;
    result.alignedB = alignedB;

    
    result.score = bestScore;

    return result;
}

TrieNode* buildTrie(const vector<string>& motifs) {
    TrieNode* root = new TrieNode();

    for (const string& motif : motifs) {
        TrieNode* current = root;

        for (char c : motif) {
            // If no path exists for this character, create one
            if (current->children.count(c) == 0) {
                current->children[c] = new TrieNode();
            }
            // Follow the path
            current = current->children[c];
        }
        // Mark end of motif
        current->isEndOfWord = true;
    }

    return root;
}

vector<string> searchMotifs(const string& protein, TrieNode* root) {
    vector<string> found;

    // Try starting a search from every position in the protein
    for (size_t i = 0; i < protein.size(); i++) {
        TrieNode* current = root;
        string match = "";

        // Walk the trie from this position
        for (size_t j = i; j < protein.size(); j++) {
            char c = protein[j];

            // No path for this character — stop searching from position i
            if (current->children.count(c) == 0) break;

            // Follow the path and build the match string
            current = current->children[c];
            match += c;

            // If we're at the end of a motif, record it
            if (current->isEndOfWord) {
                found.push_back(match);
            }
        }
    }

    return found;
}

// Clean up heap memory when done
void deleteTrie(TrieNode* root) {
    for (auto& pair : root->children) {
        deleteTrie(pair.second);
    }
    delete root;
}
// Take a list of motifs, insert each one into the trie character by character
vector<string> searchMotifs(const string& protein, TrieNode* root);
// traverse through protein string, and at every starting position, follow ther trie as far as possible, recording any motifs

