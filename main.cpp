#include <iostream>
#include "fasta.h"
#include <string>
#include <fstream>
#include <mutex>
#include <thread>
#include <functional>
#include <vector>
using namespace std;



int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: ./analyzer <filename.fasta> [motif1 motif2 ...]" << endl;
        return 1;
    }

    string filename = argv[1];
    vector<string> motifs;
    for (int i = 2; i < argc; i++) {
        motifs.push_back(argv[i]);
    }

    vector<FastaRecord> records = parseFasta(filename);

    vector<thread> threads;
    for (const auto& rec : records) {
        threads.emplace_back(analyzeRecord, ref(rec), motifs);
    }
    for (auto& t : threads) {
        t.join();
    }
    return 0;
}