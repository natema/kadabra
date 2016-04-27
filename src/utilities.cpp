#include "utilities.h"
#include <random>
#include <limits.h>
#include <sys/time.h>
using namespace std;

double get_time_sec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec + tv.tv_usec / 1000000.0 ;
}
