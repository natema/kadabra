#include "Ranking_list.h"

#include<cstdlib>
#include<climits>
#include<iostream>

// Creates a new Ranking_list of length k.
Ranking_list::Ranking_list( const uint32_t k ) : k(k)
{
    elements = (uint32_t*) malloc(k * sizeof(uint32_t));
    values = (double*) malloc(k * sizeof(double));

    for (uint32_t i = 0; i < k; i++) {
        elements[i] = i;
        values[i] = 0;
    }
}

// Updates the value of element with new_value, if already present; otherwise, it inserts it.
void Ranking_list::put(const uint32_t element, const double new_value) {
    uint32_t cur_element = element, tmp_element;
    double cur_value = new_value, tmp_value;
    uint32_t i;
    for (i = k; i-- > 0; ) {
        if (values[i] > new_value) {
            break;
        }
    }
    i++;
    while (i < k) {
        tmp_element = cur_element;
        cur_element = elements[i];
        elements[i] = tmp_element;
        tmp_value = cur_value;
        cur_value = values[i];
        values[i] = tmp_value;

        if (cur_element == element) {
            break;
        }
        i++;
    }
}


Ranking_list::~Ranking_list()
{
    free(elements);
    free(values);
}
