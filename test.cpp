#include <iostream>
#include <cassert>
#include "fasta.h"
using namespace std;

// Real sequence from JX432019.1 - Homo sapiens IGHV gene
const string REAL_DNA = 
    "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCGGTGAAGGTCTCCTGCAAGG"
    "CTTCTGGAGGCACCTTCAGCAGGTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTG"
    "GATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATT"
    "ACCGCGGACGAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGT"
    "ATTACTGTGCGAGAGATGCTCGGTTCTCTTCTAGTTGGAGTGGCCGCTTCTACTACTACATGGACGTCTG"
    "G";

void testTranscribeDNA() {
    // Basic sanity check
    assert(transcribeDNA("ATCG") == "AUCG");

    // Verify no T remains in real sequence after transcription
    string rna = transcribeDNA(REAL_DNA);
    assert(rna.find('T') == string::npos);
    assert(rna.find('U') != string::npos);

    cout << "testTranscribeDNA passed" << endl;
}

void testTranslateRNA() {
    // Basic sanity check
    string rna = "AUGUUUUAA";  // MET, PHE, STOP
    assert(translateRNA(rna) == "MF");

    // Verify real sequence translates without crashing
    // and produces a non-empty protein
    string realRna = transcribeDNA(REAL_DNA);
    string protein = translateRNA(realRna);
    assert(!protein.empty());

    // First codon of real sequence is CAG → Q (GLN)
    assert(protein[0] == 'Q');

    cout << "testTranslateRNA passed" << endl;
}

void testComputeStats() {
    // Basic sanity check
    string protein = "MMA";
    ProteinStats stats = computeStats(protein);
    assert(stats.length == 3);
    assert(stats.aminoAcidComposition["M"] == 2);
    assert(stats.aminoAcidComposition["A"] == 1);

    // Test on real translated sequence
    string realRna = transcribeDNA(REAL_DNA);
    string realProtein = translateRNA(realRna);
    ProteinStats realStats = computeStats(realProtein);

    assert(realStats.length > 0);
    assert(realStats.molecularWeight > 0.0);
    assert(!realStats.aminoAcidComposition.empty());

    cout << "testComputeStats passed" << endl;
}

void testDetectSequenceType() {
    assert(detectSequenceType(REAL_DNA) == "DNA");
    assert(detectSequenceType("AUCGAUCG") == "RNA");
    assert(detectSequenceType("MFKAGW") == "PROTEIN");

    cout << "testDetectSequenceType passed" << endl;
}

void testTrie() {
    // Basic sanity check
    vector<string> motifs = {"MF", "AG", "KL"};
    TrieNode* root = buildTrie(motifs);
    vector<string> found = searchMotifs("MFKLAG", root);
    
    cout << "  Found " << found.size() << " motifs:" << endl;
    for (const string& m : found) {
        cout << "    " << m << endl;
    }
    assert(found.size() == 3);
    assert(found[0] == "MF");
    assert(found[1] == "KL");
    assert(found[2] == "AG");
    deleteTrie(root);

    // Test on real protein sequence
    string realRna = transcribeDNA(REAL_DNA);
    string realProtein = translateRNA(realRna);
    TrieNode* realRoot = buildTrie(motifs);
    vector<string> realFound = searchMotifs(realProtein, realRoot);
    deleteTrie(realRoot);
    // Just verify it runs without crashing on a real sequence
    cout << "  Motifs found in real sequence: " << realFound.size() << endl;

    cout << "testTrie passed" << endl;
}

void testSmithWaterman() {
    // Basic known alignment
    AlignmentResult result = smithWaterman("MFKAGW", "FKA");
    assert(result.score > 0);
    assert(result.alignedA == "FKA");
    assert(result.alignedB == "FKA");

    // Test on real sequence against a known IGHV motif
    string realRna = transcribeDNA(REAL_DNA);
    string realProtein = translateRNA(realRna);
    AlignmentResult realResult = smithWaterman(realProtein, "QVQL");
    assert(realResult.score > 0);
    cout << "  Best local alignment score: " << realResult.score << endl;
    cout << "  Aligned A: " << realResult.alignedA << endl;
    cout << "  Aligned B: " << realResult.alignedB << endl;

    cout << "testSmithWaterman passed" << endl;
}

int main() {
    testTranscribeDNA();
    testTranslateRNA();
    testComputeStats();
    testDetectSequenceType();
    testTrie();
    testSmithWaterman();

    cout << "\nAll tests passed!" << endl;
    return 0;
}