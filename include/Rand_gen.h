#ifndef RAND_GEN_H
#define RAND_GEN_H

#include <cstdint>

class Rand_gen
{
    public:
        Rand_gen(const uint64_t seed);
        virtual ~Rand_gen();
        uint64_t get();
        uint64_t get_max( const uint64_t i );
    protected:
    private:
        constexpr static uint64_t a = 6364136223846793005L;
        constexpr static uint64_t c = 1442695040888963407L;
        uint64_t seed;
};

#endif // RAND_GEN_H
