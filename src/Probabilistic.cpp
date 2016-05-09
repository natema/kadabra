#include <string>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <iomanip>
#include <time.h>
#include <cfloat>
#include <omp.h>

#include "utilities.h"
#include "Probabilistic.h"
#include "Sp_sampler.h"

#define SEED 42


using namespace std;

// The status class contains the data about the k most central vertices.
Status::Status(const uint32_t k) : k(k) {
    approx_top_k = (double *) malloc( k*sizeof(double));
    top_k = (uint32_t *) malloc( k*sizeof(uint32_t) );
    finished = (bool *) malloc( k*sizeof(bool) );
    bet = (double*) malloc(k * sizeof(double));
    err_l = (double*) malloc(k * sizeof(double));
    err_u = (double*) malloc(k * sizeof(double));
}

Status::~Status() {
    free(approx_top_k);
    free(top_k);
    free(finished);
    free(bet);
    free(err_l);
    free(err_u);
}

// Creates the graph for running the approximation algorithm.
// For more information see the graph class.
Probabilistic::Probabilistic( const std::string &filename, const bool directed, const double verb ): Graph( filename, directed ), verbose(verb) {
    approx = (double *) calloc( get_nn(), sizeof(double) );
    delta_l_guess = (double *) calloc( get_nn(), sizeof(double) );
    delta_u_guess = (double *) calloc( get_nn(), sizeof(double) );
    time_bfs = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_comp_finished = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_critical = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    n_pairs = 0;
    vis_edges = 0;
    if (verbose > 0) {
        print_data();
    }
}

// Decides whether the algorithm should terminate
// INPUT: a Status object describing the current status of the algorithm
bool Probabilistic::compute_finished(Status *status) const {

    double *bet = status->bet;
    double *err_l = status->err_l;
    double *err_u = status->err_u;
    bool all_finished = true;

    uint32_t i;
    for (i = 0; i < status->k-1; i++) {
        bet[i] = status->approx_top_k[i] / status->n_pairs;
        err_l[i] = compute_f( bet[i], status->n_pairs, delta_l_guess[status->top_k[i]] );
        err_u[i] = compute_g( bet[i], status->n_pairs, delta_u_guess[status->top_k[i]] );
    }
    bet[i] = status->approx_top_k[i] / status->n_pairs;
    err_l[i] = compute_f( bet[i], status->n_pairs, this->delta_l_min_guess );
    err_u[i] = compute_g( bet[i], status->n_pairs, this->delta_u_min_guess );

    if (absolute) {
        for (uint32_t i = 0; i < status->k; i++) {
            status->finished[i] = (err_l[i] < err && err_u[i] < err);
            all_finished = all_finished && status->finished[i];
        }
    } else {
        for (uint32_t i = 0; i < status->k; i++) {
            if (i == 0) {
                status->finished[i] = (bet[i]-err_l[i] > bet[i+1]+err_u[i+1]);
            } else if (i < k) {
                status->finished[i] = (bet[i-1]-err_l[i-1] > bet[i]+err_u[i]) && (bet[i]-err_l[i] > bet[i+1]+err_u[i+1]);
            } else {
                status->finished[i] = bet[k-1]-err_u[k-1] > bet[i]+err_u[i];
            }
            status->finished[i] = status->finished[i] || (err_l[i] < err && err_u[i] < err);
            all_finished = all_finished && status->finished[i];
        }
    }

    return all_finished;
}

// Computes the function f that bounds the betweenness of a vertex from below.
// For more information, see Borassi, Natale (2016).
double Probabilistic::compute_f( const double btilde, const uint64_t iter_num, const double delta_l ) const {
    double tmp = (((double) omega) / iter_num - 1./3);
    double err_chern = (log(1./delta_l)) * 1./iter_num * (-tmp + sqrt(tmp * tmp + 2 * btilde * omega / (log(1./delta_l))));
    return min(err_chern, btilde);
}

// Computes the function g that bounds the betweenness of a vertex from above.
// For more information, see Borassi, Natale (2016).
double Probabilistic::compute_g( const double btilde, const uint64_t iter_num, const double delta_u ) const {
    double tmp = (((double) omega) / iter_num + 1./3);
    double err_chern = (log(1./delta_u)) * 1./iter_num * (tmp + sqrt(tmp * tmp + 2 * btilde * omega / (log(1./delta_u))));
    return min(err_chern, 1-btilde);
}

// Outputs the current status.
// INPUT: a Status object describing the current status, and a flag "full".
// If full is true, we output more data.
void Probabilistic::print_status(Status *status, const bool full) const {
    if (full) {
        std::cout << std::setprecision(6) << endl << "Finished after " << status->n_pairs << " iterations." << endl;
    } else {
        std::cout << std::setprecision(6) << endl << "Situation after " << status->n_pairs << " iterations." << endl;
    }
    std::cout << "Edges visited: " << vis_edges << endl;
    std::cout << "Average edges visited: " << vis_edges/status->n_pairs << endl;
    std::cout << "Total time: " << get_time_sec() - start_time << endl;
    std::cout << "Time bfs: " << time_bfs[omp_get_thread_num()] << endl;
    std::cout << "Time critical: " << time_critical[omp_get_thread_num()] << endl;
    std::cout << "Time compute finished: " << time_comp_finished[omp_get_thread_num()] << endl;
    std::cout << "(Printing thread: " << omp_get_thread_num() << ")" << endl;

    if (absolute) {
        double max_interval = 0;
        for (uint32_t i = 0; i < status->k; i++) {
            uint32_t v = status->top_k[i];
            max_interval = max(max_interval, compute_f(status->approx_top_k[i] / status->n_pairs, status->n_pairs, delta_l_guess[v]));
            max_interval = max(max_interval, compute_g(status->approx_top_k[i] / status->n_pairs, status->n_pairs, delta_u_guess[v]));
        }
        cout << "Maximum confidence interval: " << max_interval;
    }
    else {
        compute_finished(status);
        uint32_t i;
        for (i = 0; i < k; i++) {
            double bet = status->approx_top_k[i] / status->n_pairs;
            if (status->finished[i]) {
                cout << std::setw(8) << to_string(i+1) << ") ";
            } else {
                cout << std::setw(8) << "? " + to_string(i+1) << ") ";
            }
            cout << std::setw(8) << status->top_k[i] << " " << bet-compute_f(bet, status->n_pairs, delta_l_guess[status->top_k[i]]) << " ";
            cout << bet << " " << bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]) << endl;
        }
        if (full) {
            double betk = status->approx_top_k[k-1] / status->n_pairs;
            double lbetk = betk - compute_f(betk, status->n_pairs, delta_l_guess[status->top_k[k-1]]);
            uint32_t pos = k+1;
            for (i = k; i < status->k; i++) {
                double bet = status->approx_top_k[i] / status->n_pairs;
                if (bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]) > lbetk) {
                    cout << std::setw(8) << to_string(pos++) << ") ";
                    cout << std::setw(8) << status->top_k[i] << " " << bet-compute_f(bet, status->n_pairs, delta_l_guess[status->top_k[i]]) << " ";
                    cout << bet << " " << bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]) << endl;
                }
            }
        } else {
            double max_upper = 0;
            for (i = k; i < status->k; i++) {
                double bet = status->approx_top_k[i] / status->n_pairs;
                max_upper = max(max_upper, bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]));
            }
            double bet = status->approx_top_k[status->k-1] / status->n_pairs;
            max_upper = max(max_upper, bet+compute_g(bet, status->n_pairs, delta_u_min_guess));

            cout << std::setw(8) << "Others" << ") <" << max_upper;
        }
    }
    cout << endl;
}

// Sample one shortest path and updates the ranking of the betweenness approximations.
void Probabilistic::one_round(Sp_sampler &sp_sampler) {
    time_bfs[omp_get_thread_num()] -= get_time_sec();
    vector<uint32_t> path = sp_sampler.random_path();
    time_bfs[omp_get_thread_num()] += get_time_sec();

    time_critical[omp_get_thread_num()] -= get_time_sec();
    #pragma omp critical
    {
        n_pairs++;
        vis_edges += sp_sampler.vis_edges;

        for(uint32_t u:path){
            approx[u]++;
            top_k->put(u, approx[u]);
        }
    }

    time_critical[omp_get_thread_num()] += get_time_sec();
}

// Fills the input variable Status in a synchronized way.
void Probabilistic::get_status (Status *status) const {
    time_critical[omp_get_thread_num()] -= get_time_sec();
    #pragma omp critical
    {
        if (status != NULL) {
            for(uint32_t i=0; i<union_sample; i++) {
                status->top_k[i] = top_k->get(i);
                status->approx_top_k[i] = approx[status->top_k[i]];
            }
            status->n_pairs = n_pairs;
        }
    }
    time_critical[omp_get_thread_num()] += get_time_sec();
}

// Compute the values of err from which the *best* deltas are computed.
// The results are stored in err_l and err_u.
void Probabilistic::compute_bet_err(Status *status, double *bet, double *err_l, double *err_u) const {

    uint32_t i;
    double max_err = sqrt(start_factor) * err / 4;

    for (i = 0; i < status->k; i++) {
        bet[i] = status->approx_top_k[i] / status->n_pairs;
    }
    if (absolute) {
        for (i = 0; i < status->k; i++) {
            err_l[i] = err;
            err_u[i] = err;
        }
    } else {
        err_u[0] = max(err, (bet[0] - bet[1]) / 2.);
        err_l[0] = 10;
        for (i = 1; i < k; i++) {
            err_l[i] = max(err, (bet[i-1]-bet[i]) / 2.);
            err_u[i] = max(err, (bet[i]-bet[i+1]) / 2.);
        }
        for (i = k; i < status->k; i++) {
            err_l[i] = 10;
            err_u[i] = max(err, bet[k-1] + (bet[k-1]-bet[k]) / 2. - bet[i]);
        }
        for (i = 0; i < k-1; i++) {
            if (bet[i] - bet[i+1] < max_err) {
                err_l[i] = err;
                err_u[i] = err;
                err_l[i+1] = err;
                err_u[i+1] = err;
            }
        }
        for (i = k+1; i < status->k; i++) {
            if (bet[k] - bet[i] < max_err) {
                err_l[k] = err;
                err_u[k] = err;
                err_l[i] = err;
                err_u[i] = err;
            }
        }
    }
}

// Compute the *best* deltas for minimizing the stopping time of the algorithm.
// The computation is based on the heuristic in the paper Borassi, Natale (2016).
void Probabilistic::compute_delta_guess() {
    double balancing_factor = 0.001;
    double a = 0, b = 1. / err / err * log(get_nn() * 4 * (1-balancing_factor) / delta), c=(a+b)/2;
    double sum;
    Status status(union_sample);
    get_status(&status);

    double *bet = (double*) malloc(status.k * sizeof(double));
    double *err_l = (double*) malloc(status.k * sizeof(double));
    double *err_u = (double*) malloc(status.k * sizeof(double));

    compute_bet_err(&status, bet, err_l, err_u);
    for (uint32_t i = 0; i < union_sample; i++) {
        uint32_t v = status.top_k[i];
        approx[v] = approx[v] / n_pairs;
    }

    while (b-a > err/10) {
        c = (b+a)/2;
        sum = 0;
        for (uint32_t i = 0; i < union_sample; i++) {
            sum += exp(-c * err_l[i] * err_l[i] / bet[i]);
            sum += exp(-c * err_u[i] * err_u[i] / bet[i]);
        }
        sum += exp(-c * err_l[union_sample-1] * err_l[union_sample-1] / bet[union_sample-1]) * (get_nn() - union_sample);
        sum += exp(-c * err_u[union_sample-1] * err_u[union_sample-1] / bet[union_sample-1]) * (get_nn() - union_sample);

        if (sum >= delta / 2 *(1-balancing_factor)) {
            a = c;
        } else {
            b = c;
        }
    }
    delta_l_min_guess = exp(-b * err_l[union_sample-1] * err_l[union_sample-1] / bet[union_sample-1]) + delta * balancing_factor / 4. / get_nn();;
    delta_u_min_guess = exp(-b * err_u[union_sample-1] * err_u[union_sample-1] / bet[union_sample-1]) + delta * balancing_factor / 4. / get_nn();;

    for (uint32_t v = 0; v < get_nn(); v++) {
        delta_l_guess[v] = delta_l_min_guess;
        delta_u_guess[v] = delta_u_min_guess;
    }
    for (uint32_t i = 0; i < union_sample; i++) {
        uint32_t v = status.top_k[i];
        delta_l_guess[v] = exp(-b * err_l[i] * err_l[i] / bet[i]) + delta * balancing_factor / 4. / get_nn();;
        delta_u_guess[v] = exp(-b * err_u[i] * err_u[i] / bet[i]) + delta * balancing_factor / 4. / get_nn();;
    }
    free(bet);
    free(err_l);
    free(err_u);
}


// Runs the algorithm.
// INPUT: k is the number of betweenness that have to be approximated (if k=0 all betweenness
// are approximated with absolute error); delta is the probabilistic guarantee; err is the
// maximum error allowed; union_sample and start_factor are parameters of the algorithm
// that are automatically chosen.
void Probabilistic::run(uint32_t k, double delta, double err, uint32_t union_sample, uint32_t start_factor) {
    this->absolute = (k == 0);
    this->err = err;
    this->delta = delta;
    this->start_factor = start_factor;
    this->omega = 0.5 / err / err * (log2(estimate_diameter()-1) + 1 + log(0.5 / delta));
    uint32_t tau = omega / start_factor; // Da sistemare


    if (union_sample == 0) {
        union_sample = min(get_nn(), (uint32_t) max( 2 * sqrt(get_ne()) / omp_get_max_threads(), k+20. ));
    }
    this->union_sample=union_sample;
    this->k=min(k, get_nn());

    last_output = get_time_sec();
    start_time = get_time_sec();
    this->top_k = new Ranking_list(union_sample);
    srand( SEED );
    uint32_t *random_seed = (uint32_t *) malloc( omp_get_max_threads()*sizeof(uint32_t) );
    for( int i=0; i < omp_get_max_threads(); i++ ){
        random_seed[i] = rand();
    }

    #pragma omp parallel
    {
        Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
        while (n_pairs <= tau) {
            one_round(sp_sampler);
            double current_time = get_time_sec();

            if (verbose > 0 && current_time - last_output > verbose) {
                #pragma omp critical
                {
                    if (current_time - last_output > verbose) {
                        last_output = current_time;
                        cout << "First visits: " << n_pairs << "/" << tau << ".\n";
                    }
                }
            }
        }
    }

    *time_bfs = 0;
    *time_critical = 0;
    compute_delta_guess();
    n_pairs = 0;
    delete(this->top_k);
    this->top_k = new Ranking_list(union_sample);
    for (uint32_t i = 0; i < get_nn(); i++) {
        approx[i] = 0;
    }

    #pragma omp parallel
    {
        Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
        Status status(union_sample);
        status.n_pairs = 0;
        bool stop = false;

        while( !stop && status.n_pairs < omega ) {
            for (uint32_t i = 0; i <= 10; i++) {
                one_round(sp_sampler);
            }
            get_status (&status);
            time_comp_finished[omp_get_thread_num()] -= get_time_sec();
            stop = compute_finished(&status);
            time_comp_finished[omp_get_thread_num()] += get_time_sec();

            double current_time = get_time_sec();

            if (verbose > 0 && current_time - last_output > verbose) {
                #pragma omp critical(print)
                {
                    if (current_time - last_output > verbose) {
                        last_output = current_time;
                        print_status(&status);
                    }
                }
            }
        }
    }
    if (verbose > 0) {
        Status status(union_sample);
        get_status(&status);
        print_status(&status, true);
    }
    free(random_seed);
    n_pairs += tau;
}

// Destructor of the class Probabilistic.
Probabilistic::~Probabilistic() {
    free(approx);
    free(delta_l_guess);
    free(delta_u_guess);
    delete(top_k);
}

