#include <vector>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <set>
#include <math.h>
#include "Graph.h"
#include <limits.h>

using namespace std;

// Creates an empty graph.
// INPUT: The number of nodes of the graph and the list of corresponding maximum degrees.
// NOTE: This constructor is used for the pred graph used in the BFS.
Graph::Graph( const uint32_t nodes_num, const uint32_t *maxdegs ) {
    this->nodes_num = nodes_num;
    edges_num = 0;
    adj = (uint32_t **) malloc( nodes_num * sizeof(uint32_t *) );
    degrees = (uint32_t *) malloc( nodes_num * sizeof(uint32_t ) );
    cc = (uint32_t *) malloc( nodes_num * sizeof(uint32_t ) );
    inc = NULL;
    in_degrees = NULL;

    for( uint32_t i=0; i < nodes_num ; i++ ) {
        degrees[i] = 0;
        adj[i] = (uint32_t *) malloc( maxdegs[i] * sizeof(uint32_t) );
    }
}

// Loads a graph from a file.
// INPUT: The path to the file where the graph is stored and a boolean
// which tells whether the graph is directed.
// NOTE: The input format consists of a list of edges (two space-separated
// integer values per line which are the two endpoints of the edge).
// Multiple edges and self-loops are ignored.
// This is the main costructor.
Graph::Graph( const string &graph_path, const bool directed ) : directed(directed)
{
    vector< set< uint32_t > > pre_adj(0);
    vector< set< uint32_t > > pre_inc(0);
    ifstream fin(graph_path);
    uint32_t u=0, v=0;
    string line;

    while( getline(fin, line) ){
        if( line[0] != '#' ){
            sscanf( line.c_str(), "%u %u", &u, &v );
            pre_adj.resize( max( (uint32_t) pre_adj.size(), 1+max( u, v ) ), set < uint32_t > () );
            pre_adj[u].insert(v);
            if (directed) {
                pre_inc.resize( max( (uint32_t) pre_inc.size(), 1+max( u, v ) ), set < uint32_t > () );
                pre_inc[v].insert(u);
            } else {
                pre_adj[v].insert(u);
            }
        }
    }
    fin.close();
    nodes_num = pre_adj.size();
    edges_num = 0;

    adj = (uint32_t **) malloc( nodes_num * sizeof(uint32_t *) );
    degrees = (uint32_t *) malloc( nodes_num * sizeof(uint32_t ) );
    if (directed) {
        inc = (uint32_t **) malloc( nodes_num * sizeof(uint32_t *) );
        in_degrees = (uint32_t *) malloc( nodes_num * sizeof(uint32_t ) );
    } else {
        inc = adj;
        in_degrees = degrees;
    }

    for( uint32_t i=0; i < nodes_num ; i++ ){
       edges_num += pre_adj[i].size();
       adj[i] = (uint32_t *) malloc( pre_adj[i].size() * sizeof(uint32_t) );
       uint32_t j=0;
       for( auto v:pre_adj[i] ){
           adj[i][j++] = v;
       }
       degrees[i] = j;
       if (directed) {
            inc[i] = (uint32_t *) malloc( pre_inc[i].size() * sizeof(uint32_t) );
            uint32_t j=0;
            for( auto v:pre_inc[i] ){
                inc[i][j++] = v;
            }
            in_degrees[i] = j;
       }
    }
    cc = (uint32_t *) malloc( nodes_num * sizeof(uint32_t ) );
}

// Prints some data about the graph.
void Graph::print_data() {
    if (directed) {
        cout << "Directed graph\n";
    } else {
        cout << "Undirected graph\n";
    }
    cout << "Number of nodes: " << get_nn() << endl;
    cout << "Number of edges: " << get_ne() << endl;
}

// Computes the (strongly) connected components of the input graph.
// The number of components is stored in the variable n_components,
// and the component of each vertex can be found in variable cc.
// The components are in reversed topological order.
// If n_components is not 0, it means that this function has already
// run, and no more action is performed.
void Graph::compute_cc() {
    if (n_components != 0) {
        return;
    }
    uint32_t u,v,w, n = get_nn(), current_index = 0;
    uint32_t *index = (uint32_t*) malloc(n * sizeof(uint32_t));
    uint32_t *pred = (uint32_t*) malloc(n * sizeof(uint32_t));
    uint32_t *lowlink = (uint32_t*) malloc(n * sizeof(uint32_t));
    uint32_t *dfs_stack = (uint32_t*) malloc((get_ne() + 1) * sizeof(uint32_t));
    uint32_t scc_stack_end, dfs_stack_end;
    uint32_t *scc_stack = (uint32_t*) malloc(n * sizeof(uint32_t));
    //Used to keep track of which nodes are in the "current" SCC
    bool *in_scc_stack = (bool *) calloc(n, sizeof(bool));
    short *visited = (short *) calloc(n, sizeof(short));
    // The variable visited[v] is 0 if the vertex has never been visited, 1 if
    // it is an ancestor of the current vertex, 2 otherwise.
    for (u = 0; u < n; u++) {
        if (visited[u] != 0) {
            continue;
        }
        // Perform a DFS from u
        dfs_stack_end = 1;
        scc_stack_end = 0;
        dfs_stack[0] = u;
        pred[u] = u;

        while (dfs_stack_end > 0) {
             v = dfs_stack[dfs_stack_end - 1];
             if (visited[v] == 0) {
                 // It means that this is the first time we visit v.
                 // We set the index and the lowlink to be equal: during the
                 // algorithm, the lowlink may decrease.
                 visited[v] = 1;
                 index[v] = current_index;
                 lowlink[v] = current_index;
                 current_index++;
                 // We add v to the stack of vertices in the current SCC
                 scc_stack[scc_stack_end++] = v;
                 in_scc_stack[v] = 1;

                 // We iterate over all neighbors of v
                 for (uint32_t i = 0; i < degrees[v]; i++) {
                     w = adj[v][i];
                     if (visited[w] == 0) {
                         // Vertex w is added to the DFS stack
                         pred[w] = v;
                         dfs_stack[dfs_stack_end++] = w;
                     }
                     else if (in_scc_stack[w]) {
                         // We update the lowlink of v (later, we will "pass"
                         // this updated value to all ancestors of v.
                         lowlink[v] = min(lowlink[v], lowlink[w]);
                     }
                 }
             }
             else {
                 // The vertex v has already been visited.
                 dfs_stack_end--;

                 if (visited[v] == 1) {
                     // It means that we have just processed all the DFS
                     // subtree rooted at v. Hence, the lowlink of v is the
                     // final value, and we "pass" this value to the
                     // predecessor of v.
                     lowlink[pred[v]] = min(lowlink[pred[v]], lowlink[v]);

                     if (lowlink[v] == index[v]) {
                         // The DFS subtree rooted at v is a new SCC. We
                         // recover the SCC from scc_stack.
                         w = -1;
                         while (w != v) {
                             scc_stack_end--;
                             w = scc_stack[scc_stack_end];
                             in_scc_stack[w] = 0;
                             cc[w] = n_components;
                         }
                         n_components++;
                     }
                     visited[v] = 2;
                 }
            }
        }
    }
    free(index);
    free(pred);
    free(lowlink);
    free(dfs_stack);
    free(scc_stack);
    free(in_scc_stack);
    free(visited);
}

// For each strongly connected component, chooses a "pivot" vertex.
// The pivot vertex is the vertex maximizing the sum of the out-degree
// and the in-degree. The result is stored in the array pivots, which
// should already be allocated and should have length n_components.
void Graph::compute_pivot( uint32_t *pivots ) {
    for( uint32_t i=0; i<n_components; i++ ) {
        pivots[i] = UINT_MAX;
    }
    for( uint32_t i=0; i<get_nn(); i++ ) {
        uint32_t *p = &(pivots[cc[i]]);
        if ( *p==UINT_MAX || (degrees[i]+in_degrees[i]) > degrees[*p]+in_degrees[*p] ) {
            *p = i;
        }
    }
}

// Computes the strongly connected component graph.
// The result is stored in the variable cc_adj.
void Graph::compute_cc_adj() {
    cc_adj.resize(n_components);

    for( uint32_t i=0; i<get_nn(); i++ ) {
        for( uint32_t j=0; j<degrees[i]; j++ ) {
            if (cc[i] != cc[adj[i][j]]) {
                cc_adj[cc[i]].insert(cc[adj[i][j]]);
            }
        }
    }
}

// Computes the eccentricity of the vertex start in its strongly connected component.
// INPUT: The vertex whose eccentricity has to be computed, a pointer to an already-allocated
// queue, a pointer to an already-allocated array for the distances of the vertices from start,
// a boolean indicating whether the forward or backward eccentricity is requested.
// NOTE: It is assumed that the strongly connected components have already been computed
// (see function compute_cc).
uint32_t Graph::compute_ecc_in_scc(const uint32_t start, uint32_t *q, int32_t *dist, const bool backward) {
    uint32_t start_q=0, end_q=0;
    uint32_t v, u;
    uint32_t neigh_num;
    uint32_t **adj;
    uint32_t *degs;

    if (backward) {
        adj = this->inc;
        degs = this->in_degrees;
    } else {
        adj = this->adj;
        degs = this->degrees;
    }

    q[end_q++]=start;
    dist[start] = 0;

    while( start_q < end_q ){
        u = q[start_q++];
        neigh_num = degs[u];

        for( uint32_t i=0; i<neigh_num; i++ ){
            v = adj[u][i];
            if( dist[v] == -1 && cc[v] == cc[u] ){
                dist[v] = dist[u] + 1;
                q[end_q++] = v;
            }
        }
    }
    return dist[q[end_q-1]];
}

// Returns an upper bound on the diameter of the graph.
// The function uses the AllCCUpperBound technique in Borassi et al. 2015.
uint32_t Graph::estimate_diameter() {
    compute_cc();
    compute_cc_adj();

    uint32_t *q = (uint32_t*) malloc(get_nn() * sizeof(uint32_t));
    int32_t *dist_f = (int32_t*) malloc(get_nn() * sizeof(int32_t));
    int32_t *dist_b = (int32_t*) malloc(get_nn() * sizeof(int32_t));
    uint32_t *pivots = (uint32_t*) malloc(n_components * sizeof(uint32_t));
    uint32_t *ecc_f_pivots_scc = (uint32_t*) malloc(n_components * sizeof(uint32_t));
    uint32_t *ecc_b_pivots_scc = (uint32_t*) malloc(n_components * sizeof(uint32_t));
    uint32_t *ecc_f_pivots = (uint32_t*) malloc(n_components * sizeof(uint32_t));
    uint32_t diam = 0;

    for (uint32_t v = 0; v < get_nn(); v++) {
        dist_f[v] = -1;
        dist_b[v] = -1;
    }

    compute_pivot(pivots);

    for (uint32_t i = 0; i < n_components; i++) {
        ecc_f_pivots_scc[i] = compute_ecc_in_scc(pivots[i], q, dist_f, false);
        ecc_b_pivots_scc[i] = compute_ecc_in_scc(pivots[i], q, dist_b, true);
    }

    for (uint32_t i = 0; i < n_components; i++) {
        ecc_f_pivots[i] = ecc_f_pivots_scc[i];
        for (auto cc : cc_adj[i]) {
            ecc_f_pivots[i] = max(ecc_f_pivots[i],
                                  ecc_f_pivots_scc[i] + 1 + ecc_b_pivots_scc[cc]
                                  + ecc_f_pivots[cc]);
        }
        diam = max(diam, ecc_f_pivots[i] + ecc_b_pivots_scc[i]);
    }

    free(q);
    free(dist_f);
    free(dist_b);
    free(pivots);
    free(ecc_f_pivots_scc);
    free(ecc_b_pivots_scc);

    return diam;
}

// Add the edges (u,v) to the adjacency list of the graph.
// Attention: this way, it is possible to add multiple edges
void Graph::add_edge( const uint32_t u, const uint32_t v ) {
    adj[u][degrees[u]++] = v;
}

// Clean the adjacency list of the graph.
// Instead of de-allocating the array adj, the function set to zero the length of each entry.
void Graph::remove_all_edges() {
    for (uint32_t v = 0; v < get_nn(); v++) {
        degrees[v] = 0;
    }
}

// Removes some edges from the adjacency list of the graph.
// Takes as input the array vertices and its length, and set to zero the degree of those
// vertices in the adjacency list.
void Graph::remove_some_edges( const uint32_t *vertices, const uint32_t length) {
    for (uint32_t i = 0; i < length; i++) {
        degrees[vertices[i]] = 0;
    }
}

// Destructor of the graph class.
Graph::~Graph()
{
    for( uint32_t i=0; i < nodes_num ; i++ ){
       free(adj[i]);
    }
    free(adj);
    free(degrees);
    free(cc);
}
