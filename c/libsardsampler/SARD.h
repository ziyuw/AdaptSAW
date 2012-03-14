#ifndef __SARD_H_
#define __SARD_H_

#include "SAWHeaps.h"
#include "GraphNodes.h"

enum MoveType {
    LL, LH, HL
};

/* Can use this as a "chain" in PT */
struct _SARDState { 
    int* S;
    double* h_effective;
    double E;
};
typedef struct _SARDState SARDState;

/* Data used for sampling and computation of MH ratios */
struct _SAWData {
    /* Depends on beta of course! */ 
    pSAWHeap p_saw_heap_L;
    pSAWHeap p_saw_heap_H;

    int* sigma_lengths; // iters_per_move X 1
    int* rho_lengths; // // iters_per_move X 1

    int** sigmas; // iters_per_move X SAWlength_max
    int** rhos; // iters_per_move X SAWlength_max

    enum MoveType move_type;

    // Mixture components of proposals
    double log_f_fwd_LL;
    double log_f_fwd_LH;
    double log_f_fwd_HL;

    double log_f_rev_LL;
    double log_f_rev_LH;
    double log_f_rev_HL;
};
typedef struct _SAWData SAWData;

struct _SARDInputModel {
    int num_nodes;
    int num_edges;
    const GraphNodes* graph_nodes;
    const int** graph_edges;
    const double* J;
    const double* h;
};
typedef struct _SARDInputModel SARDInputModel;

struct _SARDSimParams {
    double beta_true;

    double gamma_hi_max, gamma_hi_min;
    double gamma_hi;
    double gamma_lo_max, gamma_lo_min;
    double gamma_lo;

    const int* S_in;
    const int* S_reset;

    double E_in;
    int n_moves;
    int n_adaptations, nruns_per_param_update;
    int num_nodes;

    const double* h_effective_in;
    const double* h_effective_reset;

    int iters_per_move_min, iters_per_move_max;
    int iters_per_move;

    int max_SAW_length_bound_max, max_SAW_length_bound_min;
    int SAW_length_max;

    int min_SAW_length_bound_max, min_SAW_length_bound_min;
    int SAW_length_min;

    double P_LL;
    double P_LH;
    double P_HL; // should be 1-(P_LL+P_LH)
};
typedef struct _SARDSimParams SARDSimParams;

#endif // __SARD_H_

