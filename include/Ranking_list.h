#ifndef DUMB_LIST_H
#define DUMB_LIST_H

#include <cstdint>

class Ranking_list
{
    public:
        void print_list();
        double get_value( const uint32_t i ) const {return values[i];}
        Ranking_list( const uint32_t k);
        virtual ~Ranking_list();
        void put(const uint32_t element, const double new_value);
        uint32_t get(int i) { return elements[i]; }
        const uint32_t k;
    protected:
    private:
        uint32_t *elements;
        double *values;
};

#endif // DUMB_LIST_H
