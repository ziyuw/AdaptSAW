#ifndef __IM_H_
#define __IM_H_

#include "SAWHeaps.h"
#include "GraphNodes.h"
#include "IMOutData.h"

/* Can use this as a "chain" in PT */
struct _IMState {
    int* S;
    double* h_effective;
    double E;
};
typedef struct _IMState IMState;

/* Data used for sampling and computation of MH ratios */
struct _SAWData {
    /* Depends on beta of course! */ 
    pSAWHeap p_down_heap;
    pSAWHeap p_up_heap;

    int SAW_length;
    double gamma_down;
    double gamma_up;

    int* sigma;
    int* rho;

    // The "set indices" in [0...num_nodes-1]
    int n_initial;
    int n_current;

    // Proposals:
    double log_f_fwd;
    double log_f_rev;
};
typedef struct _SAWData SAWData;

struct _IMInputModel {
    int num_nodes;
    int num_edges;
    const GraphNodes* graph_nodes;
    const int** graph_edges;
    const double* J;
    const double* h;
};
typedef struct _IMInputModel IMInputModel;

struct _IMSimParams {
    double beta_true;

    const int* S_ref;
    const int* S_in;
    const int* S_reset;

    const double* h_effective_in;
    const double* h_effective_reset;

    double E_in;
    int n_moves;

    int n_adaptations;

    const double *gammas; // IMRun_Expert just uses the first gamma
    int n_gammas;

    // Defines the uppper and lower bounds of SAW length. The exact SAW length
    // value used is chosen by adaptive MCMC.
    int SAW_length_min;
    int SAW_length_max;

    int nruns_per_param_update;
};
typedef struct _IMSimParams IMSimParams;

#endif // __IM_H_

