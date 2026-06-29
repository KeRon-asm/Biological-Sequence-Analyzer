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
        cerr << "Usage: ./analyzer <filename.fasta> [query_sequence] [motif1 motif2 ...]" << endl;
        return 1;
    }

    string filename = argv[1];
    string querySequence = argc >= 3 ? argv[2] : "";
    vector<string> motifs = {"MF", "AG", "KL"};
    for (int i = 3; i < argc; i++) {
        motifs.push_back(argv[i]);
    }

    vector<FastaRecord> records = parseFasta(filename);

    vector<thread> threads;
    for (const auto& rec : records) {
        threads.emplace_back(analyzeRecord, ref(rec), ref(motifs), ref(querySequence));
    }
    for (auto& t : threads) {
        t.join();
    }

    return 0;
}


/*
# Without alignment
./analyzer EXAMPLE_1.fasta

# With alignment against a query
./analyzer EXAMPLE_1.fasta QVQL

# With alignment and custom motifs
./analyzer EXAMPLE_1.fasta QVQL MF AG WR
*/