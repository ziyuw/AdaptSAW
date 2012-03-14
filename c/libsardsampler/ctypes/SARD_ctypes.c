#include "SARDOutData.h"
#include "SARD.h"
#include "SARD_unif.h"
#include "SARD_ego.h"
#include "SARDUtils.h"
#include "GraphNodes.h"
#include "GraphNodes_ctypes.h"

#include "utils.h"

#include "ego_callbacks.h"

#include <stdio.h>

void build_input_model_ctypes
    (int num_nodes, int num_edges,
     double *J, double *h,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     SARDInputModel *in_model);

void populate_simulation_parameters_ctypes
    (double *h, int *S_in, int *S_reset,
     double *h_effective_in, double *h_effective_reset,
     double beta_true,
     double gamma_hi_max, double gamma_hi_min,
     double gamma_lo_max, double gamma_lo_min,
     int max_SAW_length_bound_max, int max_SAW_length_bound_min,
     int min_SAW_length_bound_max, int min_SAW_length_bound_min,
     int iters_per_move_min, int iters_per_move_max,
     int nmoves, int n_adaptations, int nruns_per_param_update, int num_nodes,
     SARDSimParams *sim_params);

void set_up_common_SARDRun_output_args_ctypes
    (double *E_samples_out, double *E_proposed_out,
     double *log_f_fwd_out, double *log_f_rev_out,
     double *MH_ratio_out, int *move_type_out,
     double *rewards_out,
     int *used_sigma_lengths_out,
     int *used_rho_lengths_out, int *iters_per_move_out,
     double *used_gammas_his, double *used_gammas_los,
     double *P_LL_out, double *P_LH_out, double *P_HL_out, int *S_end,
     SARDMCMCDiagnostics *out_data);

void free_input_model(SARDInputModel *in_model);

__attribute__((visibility("default")))
void sampleSARD_unif_ext
    (int num_nodes, int num_edges,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     double *J, double *h, int *S_in, int *S_end,
     double *h_effective_in,
     double beta_true,
     double gamma_hi_max, double gamma_hi_min,
     double gamma_lo_max, double gamma_lo_min,
     int max_SAW_length_bound_max, int max_SAW_length_bound_min,
     int min_SAW_length_bound_max, int min_SAW_length_bound_min,
     int iters_per_move_min, int iters_per_move_max,
     int nmoves, int nruns_per_param_update,
     return_sample_fn_t return_sample_fn,
     double P_LL, double P_HL, int strategy,
     double *E_samples_out, double *E_proposed_out,
     double *log_f_fwd_out, double *log_f_rev_out,
     double *MH_ratio_out, int *move_type_out,
     int *used_sigma_lengths_out, int *used_rho_lengths_out,
     int *iters_per_move_out,
     double *used_gammas_his, double *used_gammas_los,
     double *P_LL_out, double *P_LH_out, double *P_HL_out)
{
    SARDInputModel in_model;
    SARDSimParams sim_params;
    SARDMCMCDiagnostics out_data;

    build_input_model_ctypes
        (num_nodes, num_edges, J, h,
         nbrs_lists, nbrs_list_lens,
         incident_edge_lists, incident_edge_list_lens, edges,
         &in_model);

    populate_simulation_parameters_ctypes
        (h, S_in, NULL,
         h_effective_in, NULL,
         beta_true,
         gamma_hi_max, gamma_hi_min,
         gamma_lo_max, gamma_lo_min,
         max_SAW_length_bound_max, max_SAW_length_bound_min,
         min_SAW_length_bound_max, min_SAW_length_bound_min,
         iters_per_move_min, iters_per_move_max,
         nmoves, 0, 0, num_nodes,
         &sim_params);

    set_up_common_SARDRun_output_args_ctypes
        (E_samples_out, E_proposed_out,
         log_f_fwd_out, log_f_rev_out,
         MH_ratio_out, move_type_out,
         NULL,
         used_sigma_lengths_out, used_rho_lengths_out, iters_per_move_out,
         used_gammas_his, used_gammas_los,
         P_LL_out, P_LH_out, P_HL_out, S_end,
         &out_data);

    sampleSARD_unif(&in_model, &sim_params, &out_data, return_sample_fn, P_LL, P_HL, strategy);

    free_input_model(&in_model);
}

__attribute__((visibility("default")))
void sampleSARD_ego_ext
    (int num_nodes, int num_edges,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     double *J, double *h,
     int *S_in, int *S_reset, int *S_end,
     double *h_effective_in,
     double *h_effective_reset,
     double beta_true,
     double gamma_hi_max, double gamma_hi_min,
     double gamma_lo_max, double gamma_lo_min,
     int max_SAW_length_bound_max, int max_SAW_length_bound_min,
     int min_SAW_length_bound_max, int min_SAW_length_bound_min,
     int iters_per_move_min, int iters_per_move_max,
     int nmoves, int n_adaptations, int nruns_per_param_update,
     acq_maximizer_fn_t acq_maximizer_fn,
     gp_add_data_fn_t add_data_to_gp_fn,
     gp_maximizer_fn_t maximize_gp_fn,
     gp_posterior_mean_fn_t gp_posterior_mean_fn,
     gp_maximizer_local_fn_t gp_maximizer_local_fn,
     sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn,
     return_sample_fn_t return_sample_fn, int return_samples,
     int strategy, double *E_samples_out, double *E_proposed_out,
     double *log_f_fwd_out, double *log_f_rev_out,
     double *MH_ratio_out, int *move_type_out,
     int *used_sigma_lengths_out, int *used_rho_lengths_out,
     int *iters_per_move_out,
     double *used_gammas_his, double *used_gammas_los,
     double *P_LL_out, double *P_LH_out, double *P_HL_out,
     double *rewards, double *policy, int policy_size)
{
    SARDInputModel in_model;
    SARDSimParams sim_params;
    SARDMCMCDiagnostics out_data;

    build_input_model_ctypes
        (num_nodes, num_edges, J, h,
         nbrs_lists, nbrs_list_lens,
         incident_edge_lists, incident_edge_list_lens, edges,
         &in_model);

    populate_simulation_parameters_ctypes
        (h, S_in, S_reset,
         h_effective_in, h_effective_reset,
         beta_true,
         gamma_hi_max, gamma_hi_min,
         gamma_lo_max, gamma_lo_min,
         max_SAW_length_bound_max, max_SAW_length_bound_min,
         min_SAW_length_bound_max, min_SAW_length_bound_min,
         iters_per_move_min, iters_per_move_max,
         nmoves, n_adaptations, nruns_per_param_update, num_nodes,
         &sim_params);

    set_up_common_SARDRun_output_args_ctypes
        (E_samples_out, E_proposed_out,
         log_f_fwd_out, log_f_rev_out,
         MH_ratio_out, move_type_out,
         rewards,
         used_sigma_lengths_out, used_rho_lengths_out, iters_per_move_out,
         used_gammas_his, used_gammas_los,
         P_LL_out, P_LH_out, P_HL_out, S_end,
         &out_data);

    sampleSARD_ego(&in_model, &sim_params, &out_data,
                   acq_maximizer_fn, add_data_to_gp_fn, maximize_gp_fn,
                   gp_posterior_mean_fn, gp_maximizer_local_fn, 
		   sample_gp_posterior_mean_fn, return_sample_fn, return_samples, strategy, policy, (size_t)policy_size);

    free_input_model(&in_model);
}

void build_input_model_ctypes
    (int num_nodes, int num_edges,
     double *J, double *h,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     SARDInputModel *in_model)
{
    GraphNodes* graph_nodes =
        build_graphnodes_ctypes
            (num_nodes, nbrs_lists, nbrs_list_lens,
             incident_edge_lists, incident_edge_list_lens);

    int** graph_edges = build_graphedges_ctypes(edges, num_edges);

    in_model->num_nodes = (int)num_nodes;
    in_model->num_edges = (int)num_edges;
    in_model->graph_nodes = (const GraphNodes*)graph_nodes;
    in_model->graph_edges = (const int**)graph_edges;
    in_model->J = J;
    in_model->h = h;
}

void populate_simulation_parameters_ctypes
    (double *h, int *S_in, int *S_reset,
     double *h_effective_in, double *h_effective_reset,
     double beta_true,
     double gamma_hi_max, double gamma_hi_min,
     double gamma_lo_max, double gamma_lo_min,
     int max_SAW_length_bound_max, int max_SAW_length_bound_min,
     int min_SAW_length_bound_max, int min_SAW_length_bound_min,
     int iters_per_move_min, int iters_per_move_max,
     int nmoves, int n_adaptations, int nruns_per_param_update, int num_nodes,
     SARDSimParams *sim_params)
{
    sim_params->S_in = (const int*)S_in;
    sim_params->S_reset = (const int*)S_reset;

    sim_params->h_effective_in = (const double*)h_effective_in;
    sim_params->h_effective_reset = (const double*)h_effective_reset;

    sim_params->n_moves = nmoves;
    sim_params->n_adaptations = n_adaptations;
    sim_params->nruns_per_param_update = nruns_per_param_update;
    sim_params->num_nodes = num_nodes;

    sim_params->E_in =
        compute_E_in(sim_params->S_in, sim_params->h_effective_in,
                     h, (int)num_nodes);

    sim_params->beta_true = beta_true;

    sim_params->gamma_hi_max = gamma_hi_max;
    sim_params->gamma_hi_min = gamma_hi_min;
    sim_params->gamma_lo_max = gamma_lo_max;
    sim_params->gamma_lo_min = gamma_lo_min;

    sim_params->max_SAW_length_bound_max = max_SAW_length_bound_max;
    sim_params->max_SAW_length_bound_min = max_SAW_length_bound_min;
    sim_params->min_SAW_length_bound_max = min_SAW_length_bound_max;
    sim_params->min_SAW_length_bound_min = min_SAW_length_bound_min;

    sim_params->iters_per_move_max = iters_per_move_max;
    sim_params->iters_per_move_min = iters_per_move_min;
}

void set_up_common_SARDRun_output_args_ctypes
    (double *E_samples_out, double *E_proposed_out,
     double *log_f_fwd_out, double *log_f_rev_out,
     double *MH_ratio_out, int *move_type_out,
     double *rewards_out,
     int *used_sigma_lengths_out, int *used_rho_lengths_out,
     int *iters_per_move_out,
     double *used_gammas_his, double *used_gammas_los,
     double *P_LL_out, double *P_LH_out, double *P_HL_out, int *S_end,
     SARDMCMCDiagnostics *out_data)
{
    out_data->E_samples = E_samples_out;
    out_data->E_proposed = E_proposed_out;
    out_data->log_f_fwd = log_f_fwd_out;
    out_data->log_f_rev = log_f_rev_out;
    out_data->MH_ratio = MH_ratio_out;
    out_data->move_type = move_type_out;

    out_data->used_sigma_lengths = used_sigma_lengths_out;
    out_data->used_rho_lengths = used_rho_lengths_out;
    out_data->iters_per_move = iters_per_move_out;

    out_data->used_gammas_his = used_gammas_his;
    out_data->used_gammas_los = used_gammas_los;

    out_data->used_P_LL = P_LL_out;
    out_data->used_P_LH = P_LH_out;
    out_data->used_P_HL = P_HL_out;

    out_data->rewards = rewards_out;

    out_data->S_end = S_end;
}

void free_input_model(SARDInputModel *in_model)
{
    destroy_graphedges(in_model->graph_edges, in_model->num_edges);
    destroy_graphnodes(in_model->graph_nodes);

    in_model->graph_edges = NULL;
    in_model->graph_nodes = NULL;
}

