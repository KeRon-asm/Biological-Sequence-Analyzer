#include <iostream>
#include <cassert>
#include "fasta.h"
using namespace std;

void testTranscribeDNA() {
    string dna = "ATCG";
    string expected = "AUCG";
    assert(transcribeDNA(dna) == expected);
    cout << "testTranscribeDNA passed" << endl;
}

void testTranslateRNA() {
    string rna = "AUGUUUUAA";  // MET, PHE, STOP
    string expected = "MF";
    assert(translateRNA(rna) == expected);
    cout << "testTranslateRNA passed" << endl;
}

void testComputeStats() {
    string protein = "MMA";  // 2x MET, 1x ALA
    ProteinStats stats = computeStats(protein);

    assert(stats.length == 3);
    assert(stats.aminoAcidComposition["M"] == 2);
    assert(stats.aminoAcidComposition["A"] == 1);

    cout << "testComputeStats passed" << endl;
}

int main() {
    testTranscribeDNA();
    testTranslateRNA();
    testComputeStats();

    cout << "All tests passed!" << endl;
    return 0;
}