#include "../sr-index.hpp"
#include <vector>
#include <fstream>
#include <sys/resource.h>

using namespace std;

int main (const int argc, const char* argv[]) {
    if (argc < 5) {
        cerr << "Usage " << argv[0] << " <reads_fasta_file> <max_read_length> <query_legth> <query_count> <query_file>\n";
        exit(1);
    }

    struct rlimit stack_limit;
    getrlimit(RLIMIT_STACK, &stack_limit);
    stack_limit.rlim_cur = RLIM_INFINITY;
    if (setrlimit(RLIMIT_STACK, &stack_limit) != 0) {
        cerr << "Unlimiting stack size failed, might crash.\n";
    }
    
    string reads_fasta_file = argv[1];
    int rlen = stoi(argv[2]);
    int qlen = stoi(argv[3]);
    int qcount = stoi(argv[4]);
    string query_file = argv[5];
    SR_index index(qlen, rlen);
    index.construct(reads_fasta_file);
    
    vector <string> queries(qcount);
    ifstream query_in(query_file, ifstream::in);
    for (int i = 0; i < qcount; i++) {
        query_in >> queries [i];
    }
    
    cerr << "Queries loaded\n";
    chrono::time_point<std::chrono::system_clock> tbegin, tend;
    chrono::duration<double> elapsed;
    tbegin = chrono::system_clock::now();
    for (int i = 0; i < qcount; i++) {
        index.find_reads(queries [i], false);
    }
    tend = chrono::system_clock::now();
    elapsed = tend - tbegin;
    cout << "Querying took " << elapsed.count() << "s\n";
}
