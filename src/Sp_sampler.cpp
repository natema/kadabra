#include "Sp_sampler.h"
#include "Rand_gen.h"
#include <iostream>

#define UNVISITED 0
#define VISITED_U 1
#define VISITED_V 2

using namespace std;

// Instantiates an object, that samples from the graph g.
// INPUT: g (a graphs), seed (the seed for the random sampler).
Sp_sampler::Sp_sampler( const Graph *g, const uint32_t seed ) {
    uint32_t n = g->get_nn();
    q = (uint32_t*) malloc( n*sizeof(uint32_t));
    ball_indicator = (uint32_t*) calloc( n, sizeof(uint32_t));
    dist = (uint32_t*) malloc( n*sizeof(uint32_t));
    uint32_t *max_deg = (uint32_t*) malloc( n*sizeof(uint32_t));

    for( uint32_t i=0; i<n; i++ ){
        max_deg[i] = max( g->degrees[i], g->in_degrees[i] );
    }
    pred = new Graph( n, max_deg);
    free(max_deg);

    randgen = new Rand_gen( seed );
    n_paths = (uint64_t*) malloc( n*sizeof(uint64_t));
    this->g = g;
}

// Returns a random node.
inline uint32_t Sp_sampler::random_node() const {
    return randgen->get_max(g->get_nn());
}

// Return a random path in the graph.
// The type returned is a vector, containing the list of all vertices
// in the path. The list of these vertices is not necessarily ordered.
// We do not return the vector by reference, because the compiler optimizations
// avoid to copy the whole vector before returning it.
vector<uint32_t> Sp_sampler::random_path() {
    // Sample sp
    uint32_t end_q = 0;
    uint64_t tot_weight = 0, cur_edge = 0;
    uint32_t random_edge;
    uint32_t v = random_node(), u = random_node();
    uint32_t *degrees = g->degrees;
    uint32_t **adj = g->adj;
    vis_edges = 0;
    // Bidirect_bfs
    vector<pair< uint32_t, uint32_t > > sp_edges;
    uint32_t x, y;
    bool have_to_stop = false;
    uint32_t start_u = 0, start_v = 1, end_u = 1, end_v = 2, start_cur, end_cur, *new_end_cur;
    uint32_t sum_degs_u = 0, sum_degs_v = 0, *sum_degs_cur;
    uint32_t neigh_num;

    while (u == v) {
        v = random_node();
    }

    end_q = 2;

    q[0] = u;
    q[1] = v;

    ball_indicator[u] = VISITED_U;
    ball_indicator[v] = VISITED_V;

    n_paths[u] = 1;
    dist[u] = 0;
    n_paths[v] = 1;
    dist[v] = 0;

    while( !have_to_stop ) {
        // Decide which ball should be expanded
        if (sum_degs_u <= sum_degs_v) {
            // Expand from u
            start_cur = start_u;
            end_cur = end_u;
            start_u = end_q;
            new_end_cur = &end_u;
            end_u = end_q;
            sum_degs_u = 0;
            sum_degs_cur = &sum_degs_u;
            degrees = g->degrees;
            adj = g->adj;
        } else {
            // Expand from v
            start_cur = start_v;
            end_cur = end_v;
            start_v = end_q;
            new_end_cur = &end_v;
            end_v = end_q;
            sum_degs_v = 0;
            sum_degs_cur = &sum_degs_v;
            degrees = g->in_degrees;
            adj = g->inc;
        }

        // Expand first ball.
        while( start_cur < end_cur ) {
            x = q[start_cur++];
            neigh_num = degrees[x];

            for( uint32_t j=0; j<neigh_num; j++ ){
                vis_edges++;
                y = adj[x][j];
                if( ball_indicator[y] == 0 ) {
                    (*sum_degs_cur) += degrees[y];
                    n_paths[y] = n_paths[x];
                    ball_indicator[y] = ball_indicator[x];
                    q[end_q++] = y;
                    (*new_end_cur)++;
                    pred->add_edge(y, x);
                    dist[y] = dist[x] + 1;
                } else if (ball_indicator[y]!=ball_indicator[x]) {
                    have_to_stop = true;
                    sp_edges.push_back( pair< uint32_t, uint32_t > (x, y) );
                }
                else if (dist[y] == dist[x] + 1) {
                    n_paths[y] += n_paths[x];
                    pred->add_edge(y, x);
                }
            }
        }
        if (*sum_degs_cur == 0) {
            have_to_stop = true;
        }
    }

    if (sp_edges.size() == 0) {
        pred->remove_some_edges(q, end_q);
        // Have to reset ball_indicator!
        for (uint32_t i = 0; i < end_q; i++) {
            ball_indicator[q[i]] = UNVISITED;
        }
        return vector<uint32_t>();
    }

    for (pair<uint32_t, uint32_t> p : sp_edges) {
        tot_weight += n_paths[p.first]*n_paths[p.second];
    }

    random_edge = randgen->get_max(tot_weight);
    vector<uint32_t> path;

    for (pair<uint32_t, uint32_t> p : sp_edges) {
        cur_edge += n_paths[p.first]*n_paths[p.second];
        if (cur_edge > random_edge) {
            backtrack_path( u, v, p.first, path );
            backtrack_path( u, v, p.second, path );
            break;
        }
    }

    // Have to reset ball_indicator!
    for (uint32_t i = 0; i < end_q; i++) {
        ball_indicator[q[i]] = UNVISITED;
    }
    // Have to reset pred!
    pred->remove_some_edges(q, end_q);
    return path;
}

// Backtracks the pred graph to extract a random shortest path.
void Sp_sampler::backtrack_path( const uint32_t u, const uint32_t v, const uint32_t start, vector<uint32_t> &path ) {
    uint64_t tot_weight = n_paths[start];
    uint32_t random_pred, cur_pred = 0;
    uint32_t w = 0;

    if (start == u || start == v) {
        return;
    }

    path.push_back(start);
    random_pred = randgen->get_max(tot_weight);


    for ( uint32_t t=0; t<pred->get_deg(start); t++ ) {
        w = pred->get_adj(start)[t];
        cur_pred += n_paths[v];
        if (cur_pred > random_pred) {
            break;
        }
    }

    if( w!=u && w!=v ) {
        backtrack_path( u, v, w, path );
    }
}

Sp_sampler::~Sp_sampler()
{
    free(q);
    free(ball_indicator);
    free(dist);
    free(n_paths);
    delete(randgen);
    delete(pred);
}
