#include "graph.h"
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stack>

#ifndef RETRIES
#define RETRIES 100
#endif

#ifndef VERTICES_THRESHOLD
#define VERTICES_THRESHOLD 200
#endif

using namespace std;
//using namespace sdsl;


map<char, int_t> base_to_int = {{'A',0}, {'C',1}, {'G',2}, {'T',3}};


vector <int_t> encode (string& s, int k) {
    vector <int_t> result;
    if ((int) s.size() < k) return result;
    int_t roll = 0;
    for (int i = 0; i < k - 1; i++) {
        roll <<= 2;
        roll += base_to_int [s [i]];
    }
    for (int i = k - 1; i < (int) s.size(); i++) {
        roll <<= 2;
        roll += base_to_int [s [i]];
        roll &= (1LL << (2*k)) - 1;
        result.push_back(roll);
    }
    return result;
}


void Graph::add_edge(int_t from, int_t to, vector <map <int_t, int> >& edges, map<int_t, int_t>& label_compress) {
    if (label_compress.count(from) == 0) {

        // Add vertex
        label_compress [from] = number_of_vertices;
        label_decompress.push_back(from);
        edges.push_back(map <int_t, int>());
        number_of_vertices ++;
    }
    if (label_compress.count(to) == 0) {

        // Add vertex
        label_compress [to] = number_of_vertices;
        label_decompress.push_back(to);
        edges.push_back(map <int_t, int>());
        number_of_vertices ++;
    }

    // Add edge
    edges [label_compress [from]] [label_compress [to]] ++;
}

void Graph::load_edges(const string& fasta_file) {
    vector <map <int_t, int> > edges;
    map <int_t, int_t> label_compress;
    cerr << "opening " + fasta_file << endl;
    ifstream file_in(fasta_file, ifstream::in);
    string seq;
    int_t n = 0;
    while(file_in >> seq) {
        n ++;
        file_in >> seq;
        vector <int_t> encoded = encode(seq , k - 1);
        for (int j = 0; j < (int) encoded.size() - 1; j++) {
            add_edge(encoded [j], encoded [j + 1], edges, label_compress);
        }
    }
    cerr << "read " << n << " reads" << endl;
    file_in.close();
    this -> edges_for_euler.clear();
    this -> edges_for_euler.resize(number_of_vertices);
    for (int_t i = 0; i < number_of_vertices; i++) {
        for (auto v : edges [i]) {
            this -> edges_for_euler[i].push_back(this->primitive_edges);
            this -> edge_dest.push_back(v.first);
            this -> edge_count.push_back(v.second);
            this -> primitive_edges ++;
        }
        edges_for_euler [i].shrink_to_fit();
    }
    edges_for_euler.shrink_to_fit();
    cerr << "Total number of edges: " << this->primitive_edges << endl;
}

/**
 * Count the length of oriented edge from v1 to v2.
 */
int_t Graph::count_distance(int_t v1, int_t v2) {
    int_t label1 = label_decompress [v1], label2 = label_decompress [v2];
    int_t scope = (1LL << (2*k - 4)) - 1;
    label1 &= scope;
    label2 >>= 2;

    for (int i = 1; i < k+3; i++) {
        if (label1 == label2) return i;
        scope >>= 2;
        label1 &= scope;
        label2 >>= 2;
    }
    cerr << "No overlap reached\n";
    exit(1);
}

/**
 * Count the price of adding edges from first "false" vertex to
 * first "true", second to second, etc. etc.
 */
int_t Graph::count_score(vector <pair <int_t, bool> >& assignment) {
    int_t score = 0;
    unsigned int false_ = 0, true_ = 0;
    while (false_ < assignment.size()) {
        while (false_ < assignment.size() && assignment [false_].second) false_ ++;
        while (true_ < assignment.size() && !assignment [true_].second) true_ ++;
        if (false_ < assignment.size()) {
            score += this -> count_distance (assignment [false_].first, assignment [true_].first);
            false_ ++;
            true_ ++;
        }
    }
    return score;
}

long long Graph::find_edge(int from, int to) {
    for (int e : edges_for_euler [from]) {
        if (edge_dest [e] == to) return e;
    }
    cerr << "edge not found\n";
    exit(1);
}

void Graph::remove_edge(int from, long long edge_number) {
    for (int i = 0; i < edges_for_euler [from].size(); i++) {
        if (edges_for_euler [from] [i] == edge_number) {
            swap(edges_for_euler [from] [i], edges_for_euler [from].back());
            edges_for_euler [from].pop_back();
            return;
        }
    }
    cerr << "edge not found\n";
    exit(1);
}

void Graph::adjoin_edges(const string& fasta_file) {
    cerr << "adjoining edges\n";
    vector <vector <long long> > reverse_edges(this -> number_of_vertices);
    for (int i = 0; i < this -> number_of_vertices; i++) {
        for (long long e : edges_for_euler [i]) {
            reverse_edges [edge_dest [e]].push_back(e);
        }
    }

    ifstream file_in(fasta_file, ifstream::in);
    string seq;
    vector <vector <pair <long long, int> > > edge_graph(primitive_edges);

    int internal_vertex_count = 0;
    map <int_t, int> mapper;
    while (file_in >> seq) {
        file_in >> seq;
        vector <int_t> encoded = encode(seq, this -> k - 1);
        if (mapper.count(encoded [0]) == 0) {
            mapper [encoded [0]] = internal_vertex_count;
            internal_vertex_count++;
        }
        if (mapper.count(encoded [1]) == 0) {
            mapper [encoded [1]] = internal_vertex_count;
            internal_vertex_count++;
        }
        int oldvertex = mapper [encoded [1]], oldedge = find_edge(mapper [encoded [0]], mapper [encoded [1]]);
        for (int i = 2; i < (int) encoded.size(); i++) {
            if (mapper.count(encoded [i]) == 0) {
                mapper [encoded [i]] = internal_vertex_count;
                internal_vertex_count ++;
            }
            int newedge = find_edge(oldvertex, mapper [encoded [i]]), edgepos = 0;
            while (edgepos < edge_graph [oldedge].size() && edge_graph [oldedge][edgepos].first != newedge) edgepos ++;
            if (edgepos == edge_graph [oldedge].size()) edge_graph[oldedge].push_back(make_pair(newedge, 1));
            else edge_graph[oldedge][edgepos].second ++;
            oldvertex = mapper [encoded [i]];
            oldedge = newedge;
        }
    }
    file_in.close();
    cerr << "edge_graph created" << endl;
    int sumscore = 0;
    vector <long long> next_edge(primitive_edges, -1);
    vector <bool> has_previous(primitive_edges, false);

    next_edge.shrink_to_fit();
    for (int i = 0; i < number_of_vertices; i++) {
        vector <long long> left_edges = reverse_edges [i];
        vector <long long> right_edges = edges_for_euler [i];
        vector <long long> best_edges;
        int bestscore = -1;
        if (left_edges.size() >= right_edges.size()) {
            sort(left_edges.begin(), left_edges.end());
            
            do {
                int curscore = 0;
                for (int j = 0; j < (int) right_edges.size(); j++) {
                    for (int l = 0; l < edge_graph [left_edges [j]].size(); l++) {
                        if (edge_graph [left_edges [j]] [l].first == right_edges [j]) {
                            curscore += edge_graph [left_edges [j]] [l].second;
                            break;
                        }
                    }
                }
                if (curscore > bestscore) {
                    best_edges = left_edges;
                    bestscore = curscore;
                }
            } while (next_permutation(left_edges.begin(), left_edges.end()));

            for (int j = 0; j < (int) right_edges.size(); j++) {
                for (int l = 0; l < edge_graph [best_edges [j]].size(); l++) {
                    if (edge_graph [best_edges [j]] [l].first == right_edges [j]) {
                        next_edge [best_edges [j]] = right_edges [j];
                        has_previous [right_edges [j]] = true;
                        sumscore += edge_graph [best_edges [j]] [l].second;
                        break;
                    }
                }
            }
        }
        else {
            sort(right_edges.begin(), right_edges.end());

            do {
                int curscore = 0;
                for (int j = 0; j < (int) left_edges.size(); j++) {
                    for (int l = 0; l < edge_graph [left_edges [j]].size(); l++) {
                        if (edge_graph [left_edges [j]] [l].first == right_edges [j]) {
                            curscore += edge_graph [left_edges [j]] [l].second;
                            break;
                        }
                    }
                }
                if (curscore > bestscore) {
                    best_edges = right_edges;
                    bestscore = curscore;
                }
            } while (next_permutation(right_edges.begin(), right_edges.end()));
            for (int j = 0; j < (int) left_edges.size(); j++) {
                for (int l = 0; l < edge_graph [left_edges [j]].size(); l++) {
                    if (edge_graph [left_edges [j]] [l].first == best_edges [j]) {
                        next_edge [left_edges [j]] = best_edges [j];
                        has_previous [best_edges [j]] = true;
                        sumscore += edge_graph [left_edges [j]] [l].second;
                        break;
                    }
                }
            }
        }
    }

    cerr << "edges concatenated\n";
    int composite_edges = 1;
    vector <int> reverse_dest(primitive_edges);
    for (int i = 0; i < number_of_vertices; i++) {
        for (int e : edges_for_euler [i]) {
            reverse_dest [e] = i;
        }
    }
    reverse_dest.shrink_to_fit();

    for (int i = 0; i < primitive_edges; i++) {
        if (next_edge [i] >= 0 && !has_previous[i]) {
            contained_edges.push_back(vector<long long>());
            int edge = i;
            while (edge != -1) {
                remove_edge(reverse_dest [edge], edge);
                contained_edges.back().push_back(edge);
                edge = next_edge[edge];
                if (edge == i) edge = -1;
            }
            contained_edges.back().shrink_to_fit();
            edges_for_euler [reverse_dest [i]].push_back(-composite_edges);
            composite_edges ++;
        }
    }
    
}

/**
 * Tries RETRIES random assignments and picks the best one.
 * Adds picked edges to the graph, to make it eulerian.
 */
void Graph::random_assignment(vector <pair <int_t, bool> >& bad_vertices) {
    cerr << "begin asignment" << endl;
    if (bad_vertices.empty()) return;
    vector <pair <int_t, bool> > best = bad_vertices;
    int_t score = this -> count_score (best);
    srand (time(NULL) + getpid());

    for (int i = 0; i < RETRIES; i++) {
        random_shuffle (bad_vertices.begin(), bad_vertices.end());
        int_t newscore = this -> count_score (bad_vertices);
        if (newscore < score) {
            score = newscore;
            best = bad_vertices;
        }
    }
    
    unsigned int false_ = 0, true_ = 0;
    while (false_ < best.size()) {
        while (false_ < best.size() && best [false_].second) false_ ++;
        while (true_ < best.size() && !best [true_].second) true_ ++;
        if (false_ < best.size()) {
            this -> edges_for_euler [best [false_].first].push_back(this-> primitive_edges);
            this -> edges_for_euler [best [false_].first].shrink_to_fit();
            edge_dest.push_back(best [true_].first);
            edge_count.push_back(0);
            this->primitive_edges ++;
            false_ ++;
            true_ ++;
        }
    }
}

int Graph::edge_destination(long long edge) {
    if (edge >= 0) return edge_dest [edge];
    else return edge_dest [contained_edges [-1 -edge].back()];
}

void Graph::connect_components() {
    cerr << "begin connecting" << endl;

    // Some necessary structures
    vector <bool> visited(number_of_vertices, false);
    srand (time(NULL) + getpid());
    vector <vector <int> > reverse_edges(number_of_vertices);

    // Construct reverse edges for the dfs
    for (int_t i = 0; i < number_of_vertices; i++) {
        for(int edge : edges_for_euler [i]) {
            reverse_edges [edge_destination (edge)].push_back(i);
        }
    }

    int_t last_begin = -1;
    int total_connections = 0;
    for (int_t i = 0; i < number_of_vertices; i++) {
        if (!visited [i] && (reverse_edges[i].size() > 0 || edges_for_euler [i].size() > 0)) {
            if (nonempty_vertex == -1) nonempty_vertex = i;
            stack<int> buffer;
            visited [i] = true;
            buffer.push(i);
            
            // DFS
            while (!buffer.empty()) {
                int v = buffer.top();
                buffer.pop();

                // Normal edges
                for (int edge : edges_for_euler [v]) {
                    if (edge >= 0) {
                        if (!visited [edge_destination (edge)]) {
                            visited [edge_destination (edge)] = true;
                            buffer.push(edge_destination (edge));
                        }
                    }
                }

                // Reverse edges
                for (int new_v : reverse_edges [v]) {
                    if (!visited [new_v]) {
                        visited [new_v] = true;
                        buffer.push(new_v);
                    }
                }
            }
            
            // Connect
            if (last_begin != -1ull) {
                edges_for_euler [last_begin].push_back(this->primitive_edges);
                edges_for_euler [last_begin].shrink_to_fit();
                edge_count.push_back(0);
                edge_dest.push_back(i);
                this -> primitive_edges ++;
                total_connections ++;
            }
            last_begin = i;
        }
    }
    cerr << "Connected " << total_connections + 1 << "components\n";
}

void Graph::construct_edges_for_euler() {


    // Count bad vertices (outdegree != indegree);
    // second value is true iff outdegree > indegree
    // In other terms, second value is true iff we need to add inbound edges to the vertex
    vector <pair <int_t, bool> > bad_vertices;
    vector <int> degreecount (number_of_vertices, 0);
    for (int_t i = 0; i < number_of_vertices; i++) {
        degreecount [i] += this -> edges_for_euler [i].size();
        for (int edge : this -> edges_for_euler [i]) {
            degreecount [edge_destination (edge)] --;
        }
    }

    for (int_t i = 0; i < number_of_vertices; i++) {
        if (degreecount[i] > 0) {
            for(int j = 0; j < degreecount [i]; j++) {
                bad_vertices.push_back(make_pair(i, true));
            }
        }
        else {
            for (int j = degreecount [i]; j < 0; j++) {
                bad_vertices.push_back(make_pair(i, false));
            }
        }
    }
    
    cerr << "assignment problem size: " << bad_vertices.size() << endl;
    // Add edges so that the graph is eulerian.
    this->random_assignment(bad_vertices);
    // Connect all the components (in unoriented sense)
    this -> connect_components();
}

vector <int> Graph::path_counts() {
    vector <int> result_counts(1, -1000);
    for (int edge : this -> result_edges) {
        for (int primitive_edge : list_edges(edge)) {
            result_counts.push_back(edge_count [primitive_edge]);
        }
    }
    result_counts.shrink_to_fit();
    return result_counts;
}

void Graph::euler_recursive(int_t v, vector <long long>& result, int previous_edge) {
    while (!edges_for_euler[v].empty()) {
        int next = edges_for_euler [v].back();
        edges_for_euler[v].pop_back();
        euler_recursive(edge_destination (next), result, next);
    }
    result.push_back(previous_edge);
}

vector <long long> Graph::list_edges(long long edge) {
    if (edge >= 0) return vector <long long> (1, edge);
    else return contained_edges [-1 -edge];
}

vector <int_t> Graph::euler_path() {
    this -> construct_edges_for_euler();

    cerr << "begin recursive euler" << endl;


    euler_recursive(nonempty_vertex, this->result_edges, -1);
    result_edges.pop_back();
    reverse(this -> result_edges.begin(), this -> result_edges.end());
    this->result_edges.shrink_to_fit();
    vector <int_t> decompressed_result;

    decompressed_result.push_back(label_decompress[nonempty_vertex]);
    for (int edge : this -> result_edges) {
        for (int primitive_edge : list_edges(edge)) {
            decompressed_result.push_back(label_decompress [edge_dest [primitive_edge]]);
        }
    }
    decompressed_result.shrink_to_fit();
    cerr << "finish\n";
    return decompressed_result;
}

