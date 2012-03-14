#include<stdlib.h>
#include<stdbool.h>
#include "GraphNodes.h"
#include "SARDOutData.h"
#include "SAWHeaps.h"
#include "SARD.h"


void allocMCMCDiagnostics(SARDMCMCDiagnostics* p_data, int n_moves)
{
    p_data->n_moves = n_moves;

    p_data->E_samples = (double*)malloc(n_moves*sizeof(double));
    p_data->E_proposed = (double*)malloc(n_moves*sizeof(double));
    p_data->MH_ratio = (double*)malloc(n_moves*sizeof(double));
    p_data->log_f_fwd = (double*)malloc(n_moves*sizeof(double));
    p_data->log_f_rev = (double*)malloc(n_moves*sizeof(double));
    p_data->move_type = (int*)malloc(n_moves*sizeof(int));
}

void trashMCMCDiagnostics(SARDMCMCDiagnostics* p_data)
{
    free(p_data->E_samples);
    free(p_data->E_proposed);
    free(p_data->MH_ratio);
    free(p_data->log_f_fwd);
    free(p_data->log_f_rev);
    free(p_data->move_type);
}

void storeDiagnosticsData
    (SARDMCMCDiagnostics* p_data,
     const SARDSimParams *p_sim_params,
     int curr_move, double E_0, double E_prime,
     double log_f_iter_fwd, double log_f_iter_rev,
     double MH_ratio, int move_type,
     int *sigma_lengths, int *rho_lengths, int iters_per_move,
     double used_gamma_hi, double used_gamma_lo,
     double P_LL, double P_LH, double P_HL, double reward)
{
    // Note these implement COLUMN-MAJOR storage!
    const int n_moves = p_data->n_moves;

    // Generically as Matrix[curr_move][curr_replica]
    
    p_data->E_samples[curr_move] = E_0;
    p_data->E_proposed[curr_move] = E_prime;
    p_data->MH_ratio[curr_move] = MH_ratio;
    p_data->log_f_fwd[curr_move] = log_f_iter_fwd;
    p_data->log_f_rev[curr_move] = log_f_iter_rev;
    p_data->move_type[curr_move] = move_type;    
    p_data->iters_per_move[curr_move] = iters_per_move;
    
    for (int i = 0; i < iters_per_move; ++i)
    {
        int out_buf_idx = curr_move*p_sim_params->iters_per_move_max + i;
        p_data->used_sigma_lengths[out_buf_idx] = sigma_lengths[i];
        p_data->used_rho_lengths[out_buf_idx] = rho_lengths[i];
    }

    p_data->used_gammas_his[curr_move] = used_gamma_hi;
    p_data->used_gammas_los[curr_move] = used_gamma_lo;

    p_data->used_P_LL[curr_move] = P_LL;
    p_data->used_P_LH[curr_move] = P_LH;
    p_data->used_P_HL[curr_move] = P_HL;

    
    // There may not necessarily be rewards involved, i.e. non-adaptive MCMC.
    if (p_data->rewards != NULL)
    {
        p_data->rewards[curr_move] = reward;
    }
}

void storeEndState(SARDMCMCDiagnostics* p_data, int *S_end, int num_nodes)
{
    for (int i = 0; i < num_nodes; ++i)
    {
        p_data->S_end[i] = S_end[i];
    }
}

