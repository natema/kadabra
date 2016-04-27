#include <iostream>
#include <dirent.h>
#include <limits.h>
#include <math.h>
#include <ctime>
#include <cstddef>
#include <unistd.h>
#include <ctype.h>

#include "Rand_gen.h"
#include "Graph.h"
#include "utilities.h"
#include "Probabilistic.h"

extern char *optarg;
static const std::string ERROR_HEADER = "ERROR: ";
using namespace std;

bool directed = false;
double verb = 60;
double delta;
double err;
char *graph_file;
int64_t k = 0;

/**
 * Print usage on stderr.
 */
void usage(const char *binary_name) {
    std::cerr << binary_name
        << ": compute betweenness centrality approximations for all nodes"
        << std::endl;
    std::cerr << "USAGE: " << binary_name << " [-dh] [-v verbosity] [-k k_value] epsilon delta graph"
        << std::endl;
    std::cerr << "\t-d: consider the graph as directed" << std::endl;
    std::cerr << "\t-k: compute the top-k betweenness centralities (if 0, compute all of them with absolute error) " << std::endl;
    std::cerr << "\t-h: print this help message" << std::endl;
    std::cerr << "\t-v: print additional messages (verbosity is the time in second between subsequent outputs)" << std::endl;
    std::cerr << "\terr: accuracy (0 < epsilon < 1)" << std::endl;
    std::cerr << "\tdelta: confidence (0 < delta < 1)" << std::endl;
    std::cerr << "\tgraph: graph edge list file" << std::endl;
}

/**
 * Parse command line options.
 * Return 0 if everything went well, 1 if there were errors, 2 if -h was specified.
 */
int parse_command_line(int& argc, char *argv[]) {
    int opt;
    while ((opt = getopt(argc, argv, "dhk:v:")) != -1) {
        switch (opt) {
        case 'd':
            directed = true;
            break;
        case 'h':
            return 2;
            break;
        case 'k':
            k = std::strtod(optarg, NULL);
            if (errno == ERANGE || k < 0 || k > UINT_MAX) {
                std::cerr << ERROR_HEADER
                    << "The value k should be between 0 and 2^32-1."
                    << std::endl;
                return 1;
            }
            break;
        case 'v':
            verb = std::strtod(optarg, NULL);
            if (errno == ERANGE || verb < 0) {
                std::cerr << ERROR_HEADER
                    << "The verbosity should be a positive number, or 0 to produce no output."
                    << std::endl;
                return 1;
            }
            break;
        }
    }

    if (optind != argc - 3) {
        std::cerr << ERROR_HEADER << "Wrong number of arguments" << std::endl;
        return 1;
    } else {
        err = std::strtod(argv[argc - 3], NULL);
        if (errno == ERANGE || err >= 1.0 || err <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "The error err should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        delta = std::strtod(argv[argc - 2], NULL);
        if (errno == ERANGE || delta >= 1.0 || delta <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "Delta should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        graph_file = argv[argc - 1];
    }

    return 0;
}

int main(int argc, char *argv[]){
    int correct_parse = parse_command_line(argc, argv);

    if (correct_parse != 0) {
        usage(argv[0]);
        return correct_parse!=2;
    }

    Probabilistic G( graph_file, directed, verb );
    G.run((uint32_t) k, delta, err);
    return 0;
}
