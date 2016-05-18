#include "sr-index.hpp"

using namespace std;
using namespace sdsl;

map<int_t, char> int_to_base = {{0,'A'}, {1,'C'}, {2,'G'}, {3,'T'}};

string decode_two(int_t v1, int_t v2, int k, int count, vector <int>& string_counts) {
    int_t cv1 = v1, cv2 = v2;
    int_t add = 1, scope = (1LL << (2*k - 4)) - 1;
    cv1 &= scope;
    cv2 >>= 2;
    vector <int> magic;
    while (cv1 != cv2) {
        add ++;
        cv2 >>= 2;
        scope >>= 2;
        cv1 &= scope;
    }

    cv2 = v2;
    string result;
    for (int i = 0; i < add; i++) {
        result.push_back (int_to_base [cv2 & 3]);
        cv2 >>= 2;
        string_counts.push_back(count + 1);
    }
    reverse(result.begin(), result.end());
    return result;
}

string decode(vector <int_t> superstring, int k, vector <int> path_counts, vector <int>& string_counts) {
    string result;
    string_counts.clear();
    int_t label1 = superstring [0], label2;
    for (int i = 0; i < k - 1; i++) {
        result.push_back(int_to_base [label1 & 3]);
        string_counts.push_back(2);
        label1 >>= 2;
    }
    reverse(result.begin(), result.end());

    for (int i = 0; i < (int) superstring.size() - 1; i++) {
        result += decode_two(superstring [i], superstring [i+1], k, path_counts[i+1], string_counts);
    }
    return result;
}

void SR_index::construct_superstring(const string& fasta_file) {
    Graph g(this -> k);
    g.load_edges(fasta_file);
    g.adjoin_edges(fasta_file);

    vector <int_t> result_ints = g.euler_path();
    vector <int> result_counts = g.path_counts();
    vector <int> string_counts;
    
    cerr << "Construct FM-index\n";
    construct_im(this -> fm_index, decode(result_ints, this -> k, result_counts, string_counts), 1);
    
    this -> counts = vlc_vector<>(string_counts);
    cerr << fm_index.size() << ' ' << counts.size() << endl;
    cerr << "FM index size in mb: " << size_in_mega_bytes(fm_index) << endl;
    cerr << "counts vector size in mb: " << size_in_mega_bytes(counts) << endl;
    bit_vector valid_ends(string_counts.size());
    for (auto i = 0; i < string_counts.size(); i++) {
        valid_ends [i] = (string_counts [i] > 1);
    }
    valid_end = rrr_vector<>(valid_ends);
    cerr << "Valid ends vector size in mb: " << size_in_mega_bytes(valid_end) << endl;
}

void SR_index::construct(const string& fasta_file) {
    construct_superstring(fasta_file);
    cerr << "start processing intervals" << endl;
    
    ifstream file_in(fasta_file, ifstream::in);
    string seq;
    vector <long long> start_indices;
    vector <bool> starts, valid_positions;
    vector <long long> set_bits;
    long long current_offset = 0;
    int k = this -> k;
    while (file_in >> seq) {
        file_in >> seq;
        vector <pair <int, int> > positions;
        for (int i = 0; i <= (int) seq.size() - k; i++) {
            auto occs = locate(fm_index, seq.substr(i, k));
            bool found = false;
            for (int pos : occs) {
                if (counts [pos + k - 1] > 1) {
                    positions.push_back(make_pair(pos - i, i));
                    found = true;
                    break;
                }
            }
            if (!found) {
                cerr << "error" << ' ' << seq.substr(i, k) << endl;
                for (int pos : occs) {
                    cerr << pos << endl;
                }
                exit(1);
            }
        }
        vector <bool> validinread(max_read_length + 1, false); 
        sort (positions.begin(), positions.end());
        starts.push_back(true);
        start_indices.push_back(positions[0].first);
        for (int i = positions [0].second; i < positions [0].second + k; i++) validinread [i] = true;
        for (int i = 1; i < positions.size(); i++) {
            if (start_indices.back() != positions [i].first) {
                start_indices.push_back(positions [i].first);
                starts.push_back(false);
                bool last = false;
                for (int j = 0; j <= max_read_length; j++) {
                    if (validinread [j] != last) set_bits.push_back(current_offset + j); 
                    last = validinread [j];
                    validinread [j] = false;
                }
                current_offset += max_read_length + 1;
            }
            for (int j = positions [i].second; j < positions [i].second + k; j++) validinread [j] = true;
        }
        bool last = false;
        for (int j = 0; j <= max_read_length; j++) {
            if (validinread [j] != last) set_bits.push_back(current_offset + j);
            else valid_positions.push_back(false);
            last = validinread [j];
            validinread [j] = false;
        }
        current_offset += max_read_length + 1;
    }

    vector <long long> start_indices_permutation(start_indices.size());
    for (int i = 0; i < start_indices.size(); i++) {
        start_indices [i] += max_read_length;
        start_indices_permutation [i] = i;
    }
    sort (start_indices_permutation.begin(), start_indices_permutation.end(), [&start_indices](int i1, int i2) {return start_indices [i1] < start_indices [i2];});

    file_in.close();
    cerr << "start_indices number of elements: " << start_indices.size() << endl;
    this -> start_indices = vlc_vector<>(start_indices);
    cerr << "start_indices vlc: " << size_in_mega_bytes(this -> start_indices) << endl;
    this -> start_indices_permutation = vlc_vector<>(start_indices_permutation);
    cerr << "start_indices_permutation vlc " << size_in_mega_bytes(this->start_indices_permutation) << endl;
    bit_vector starts_b(starts.size()), valid_positions_b(valid_positions.size());
    for (unsigned int i = 0; i < starts.size(); i++) starts_b [i] = starts [i];
    this -> new_read_start = sd_vector<>(starts_b);
    this -> read_start_rank = sd_vector<>::rank_1_type(&(this -> new_read_start));
    cerr << "starts sd_vector: " << size_in_mega_bytes(this -> new_read_start) << endl;
    this -> valid_in_read = sd_vector<>(set_bits.begin(), set_bits.end());
    this -> valid_in_read_rank = sd_vector<>::rank_1_type(&(this -> valid_in_read));
    cerr << "valid in read sd_vector: " << size_in_mega_bytes(this -> valid_in_read) << endl;
    file_in.close();
}

vector<int> SR_index::find_reads(const string& query, bool debug) {
    if (query.size() > k) {
        cerr << "Query longer than k\n";
        exit(1);
    }
    long long the_position = -1;
    vector <int> result;
    for (auto pos : locate(fm_index, query)) {
        if (counts [pos + query.size() - 1] > 1) {
            the_position = pos;
            break;
        }
    }

    if (the_position == -1) {
        return result;
    }

    if (debug) cerr << "THE position : " << the_position << endl;

    long long lower = -1, upper = start_indices_permutation.size();
    while (upper - lower > 1) {
        int middle = (upper + lower) / 2;
        if (start_indices[start_indices_permutation[middle]]  < the_position + query.size() - 1) lower = middle;
        else upper = middle;
    }
    long long current = upper;
    if (debug) cerr << "First possible index: " << current << ' ' << start_indices [start_indices_permutation [current]] << ' ' << (long long) start_indices [start_indices_permutation [current] ] - max_read_length << endl << "Zaciatok ";
    set <int> results;
    while (current < start_indices_permutation.size() && (long long) start_indices[start_indices_permutation[current]] - max_read_length <= the_position) {
        long long curstart = start_indices [start_indices_permutation[current]] - max_read_length;
        long long valid_in_read_offset = (max_read_length + 1) * start_indices_permutation [current];
        if (debug) cerr << start_indices_permutation[current];
        if ((valid_in_read_rank(valid_in_read_offset + the_position - curstart + 1) % 2) == 1 && (valid_in_read_rank(valid_in_read_offset + the_position - curstart + query.size()) == valid_in_read_rank(valid_in_read_offset + the_position - curstart + 1))) {
            results.insert(read_start_rank(start_indices_permutation [current] + 1));
        }

        current ++;
    }
    if(debug) cerr << endl;

    for (auto x : results) {
        result.push_back(x - 1);
    }
    return result;

}

void SR_index::print_superstring() {
    cout << "superstring: " << extract(fm_index, 0, fm_index.size() - 1) << endl;
}
