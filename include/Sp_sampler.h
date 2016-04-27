#ifndef SP_SAMPLER_H
#define SP_SAMPLER_H

#include <vector>
#include <ctype.h>

#include "Rand_gen.h"
#include "Graph.h"

class Sp_sampler
{
    public:
        Sp_sampler( const Graph *g, const uint32_t seed );
        virtual ~Sp_sampler();
        std::vector<uint32_t> random_path();
        uint64_t vis_edges;
    protected:
    private:
        void backtrack_path( const uint32_t u, const uint32_t v, const uint32_t start, std::vector<uint32_t> &path );
        inline uint32_t random_node() const;
        Graph *pred;
        uint32_t *ball_indicator;
        uint32_t *dist;
        uint32_t *q;
        uint64_t *n_paths;
        Rand_gen *randgen;
        const Graph *g;
};

#endif // SP_SAMPLER_H
