#include "../sr-index.hpp"
#include <vector>
#include <unistd.h>
#include <fstream>
#include <sys/resource.h>


int main (const int argc, char* argv[]) {
    string usage = string(argv [0]) + " <k> <fasta_file>";
    if (argc < 3) {
        cout << usage << endl;
        return 1;
    }
    cerr << vector <int>().max_size() << endl;
    struct rlimit stack_limit;
    getrlimit(RLIMIT_STACK, &stack_limit);
    stack_limit.rlim_cur = RLIM_INFINITY;
    if (setrlimit(RLIMIT_STACK, &stack_limit) != 0) {
        cerr << "Unlimiting stack size failed, might crash.\n";
    }
    int k = atol(argv [1]);
    boost::filesystem::path orig_file = boost::filesystem::path(argv [2]);
    ifstream file_in(orig_file.string(), ifstream::in);
    string seq;
    file_in >> seq; file_in >> seq;
    file_in.close();
    
    SR_index index(k, seq.size());
    index.construct(orig_file.string());
    index.print_superstring();
    string query;

    file_in.open(orig_file.string(), ifstream::in);
    int counter = 0;
    while (file_in >> seq) {
        file_in >> seq;
        bool ok = false;
        string q = seq.substr(rand()%(seq.size() - k + 1), k);
        vector <int> results = index.find_reads(q, false);
        for (auto x : results) {
            if (x == counter) {
                ok = true;
                break;
            }
        }
        if (!ok) {
            cerr << "String not found where expected!\n";
            cerr << q << ' ' << counter << endl;
        }
        counter++;
    } 
    cerr << "finish" << endl;
    while (cin >> query) {
        for (auto x : index.find_reads(query, false)) {
            cout << x << ' ';
        }
        cout << endl;
    } 

} 

