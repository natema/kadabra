#include "Rand_gen.h"
#include <random>

Rand_gen::Rand_gen( uint64_t seed )
{
    this->seed = seed;
}

Rand_gen::~Rand_gen()
{
}

// Returns a random number generated according to a linear congruential generator.
// The constants are borrowed from those of MMIX by D. Knuth. 
uint64_t Rand_gen::get() {
    seed = (a * seed) + c;
    return seed;
}

// Return a random number between 0 and m-1 (included).
uint64_t Rand_gen::get_max( uint64_t m ) {
    uint32_t r;
    while( (r = get()) > (UINT64_MAX - UINT64_MAX % m) ) {}
    return r % m;
}

