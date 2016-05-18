#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <set>
#include "common.h"

using namespace std;

class Graph {
    vector <int_t> label_decompress;
    vector <vector <long long> > edges_for_euler;
    // Primitive edges
    vector <int> edge_count;
    vector <int> edge_dest;
    // Composite edges
    vector <vector <long long> > contained_edges;

    vector <long long> result_edges;
    int_t number_of_vertices = 0;
    long long primitive_edges = 0;
    int k;
    int nonempty_vertex = -1;
    void connect_components();
    void construct_edges_for_euler();
    void random_assignment(vector <pair <int_t, bool> >&);
    int_t count_distance(int_t, int_t);
    int_t count_score(vector <pair <int_t, bool> >&);
    void euler_recursive(int_t, vector <long long>&, int);
    long long find_edge(int, int);
    void remove_edge(int, long long);
    int edge_destination(long long);
    vector <long long> list_edges(long long);
    public:
        void add_edge(int_t, int_t, vector <map <int_t, int> >&, map<int_t, int_t>&);
        Graph(int kk): k(kk){}
        vector <int_t> euler_path();
        vector <int> path_counts();
        void load_edges(const string&);
        void adjoin_edges(const string&);
};

#endif //GRAPH_H
