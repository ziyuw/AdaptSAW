#include <stdlib.h>
#include "GraphNodes.h"
#include "IMOutData.h"
#include "SAWHeaps.h"
#include "IM.h"


void allocMCMCDiagnostics(IMMCMCDiagnostics* p_data, int n_moves)
{
    p_data->n_moves = n_moves;
    p_data->E_samples = (double*)malloc(n_moves*sizeof(double));
    p_data->E_proposed = (double*)malloc(n_moves*sizeof(double));
    p_data->MH_ratio = (double*)malloc(n_moves*sizeof(double));
    p_data->log_f_fwd = (double*)malloc(n_moves*sizeof(double));
    p_data->log_f_rev = (double*)malloc(n_moves*sizeof(double));

    p_data->SAW_length = (int*)malloc(n_moves * sizeof(int));
    p_data->gammas = (double*)malloc(n_moves * sizeof(double));

    p_data->rewards = (double*)malloc(n_moves*sizeof(double));
}

void trashMCMCDiagnostics(IMMCMCDiagnostics* p_data)
{
    free(p_data->E_samples);
    free(p_data->E_proposed);
    free(p_data->MH_ratio);
    free(p_data->log_f_fwd);
    free(p_data->log_f_rev);
    free(p_data->SAW_length);
    free(p_data->gammas);
    free(p_data->rewards);
}

void storeDiagnosticsData
    (IMMCMCDiagnostics* p_data,
     int curr_move, double E_0, double E_prime,
     double log_f_fwd, double log_f_rev, double MH_ratio,
     int SAW_length, double gamma)
{
    p_data->E_samples[curr_move] = E_0;
    p_data->E_proposed[curr_move] = E_prime;
    p_data->MH_ratio[curr_move] = MH_ratio;
    p_data->log_f_fwd[curr_move] = log_f_fwd;
    p_data->log_f_rev[curr_move] = log_f_rev;

    p_data->SAW_length[curr_move] = SAW_length;
    p_data->gammas[curr_move] = gamma;
}

void storeEndState(IMMCMCDiagnostics* p_data, int *S_end, int num_nodes)
{
    for (int i = 0; i < num_nodes; ++i)
    {
        p_data->S_end[i] = S_end[i];
    }
}

