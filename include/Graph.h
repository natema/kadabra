#ifndef GRAPH_H
#define GRAPH_H
#include <string>
#include <vector>
#include <set>


class Graph
{
    public:
        Graph( const uint32_t nodes_num, const uint32_t *maxdegs );
        Graph( const std::string &filename, const bool directed=false );
        inline uint32_t get_nn() const {
            return nodes_num;
        }
        const inline uint64_t get_ne(){
            return edges_num;
        }
        virtual ~Graph();
        void add_edge(const uint32_t, const uint32_t);
        void remove_all_edges();
        void remove_some_edges( const uint32_t *vertices, const uint32_t length);
        const uint32_t* get_adj( const uint32_t u) {
            return adj[u];
        }
        uint32_t get_deg(const uint32_t u) {
            return degrees[u];
        }
        uint32_t n_components=0;
        uint32_t **adj;
        uint32_t *degrees;
        uint32_t **inc;
        uint32_t *in_degrees;
        bool directed;
        uint32_t estimate_diameter();
        uint32_t *cc;
        void compute_cc();
        void print_data();
    protected:
        uint32_t nodes_num=0;
        uint64_t edges_num=0;
    private:
        void compute_pivot( uint32_t *pivots );
        std::vector < std::set < uint32_t > > cc_adj;
        uint32_t compute_ecc_in_scc(const uint32_t start, uint32_t *q, int32_t *dist, const bool backward);
        void compute_cc_adj();
};

#endif // GRAPH_H
