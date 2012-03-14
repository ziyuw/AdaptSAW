#include "IMOutData.h"
#include "IM.h"
#include "IM_unif.h"
#include "IM_ego.h"
#include "IMUtils.h"
#include "GraphNodes.h"
#include "GraphNodes_ctypes.h"

#include "utils.h"

#include "ego_callbacks.h"

void build_input_model_ctypes
    (int num_nodes, int num_edges,
     double *J, double *h,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     IMInputModel *in_model);

void populate_simulation_parameters_ctypes
    (double *h, int *S_ref, int *S_in, int *S_reset,
     double *h_effective_in, double *h_effective_reset,
     double beta_true,int len_gammas, double *gammas,
     int SAW_length_min, int SAW_length_max,
     int nmoves, int n_adaptations, int nruns_per_param_update, int num_nodes,
     IMSimParams *sim_params);

void set_up_common_IMRun_output_args_ctypes
    (double *E_samples, double *E_proposeds,
     double *log_f_fwds, double *log_f_revs,
     double *MH_ratios, int *SAW_lengths, double *used_gammas,
     double *rewards, int nactions, int *S_end,
     IMMCMCDiagnostics *out_data);

__attribute__((visibility("default")))
void sampleIM_unif_ext
    (int num_nodes, int num_edges,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     double *J, double *h, int *S_ref, int *S_in, int *S_end,
     double *h_effective_in,
     double beta_true, int len_gammas, double *gammas,
     int SAW_length_min, int SAW_length_max, int nmoves,
     int nruns_per_param_update,
     double *E_samples, double *E_proposeds,
     double *log_f_fwds, double *log_f_revs,
     double *MH_ratios, int *SAW_lengths, double *used_gammas)
{
    IMInputModel in_model;
    IMSimParams sim_params;
    IMMCMCDiagnostics out_data;

    build_input_model_ctypes
        (num_nodes, num_edges, J, h,
         nbrs_lists, nbrs_list_lens,
         incident_edge_lists, incident_edge_list_lens, edges,
         &in_model);
    populate_simulation_parameters_ctypes
        (h, S_ref, S_in, NULL, h_effective_in, NULL,
         beta_true, len_gammas, gammas,
         SAW_length_min, SAW_length_max,
         nmoves, 0, nruns_per_param_update, num_nodes,
         &sim_params);

    set_up_common_IMRun_output_args_ctypes
        (E_samples, E_proposeds, log_f_fwds, log_f_revs,
         MH_ratios, SAW_lengths, used_gammas, NULL, 0, S_end,
         &out_data);

    sampleIM_unif((const IMInputModel*)&in_model,
                   (const IMSimParams*)&sim_params,
                   &out_data);

    free_input_model(&in_model);
}
__attribute__((visibility("default")))
void sampleIM_ego_ext
    (int num_nodes, int num_edges,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     double *J, double *h,
     int *S_ref, int *S_in, int *S_reset, int *S_end,
     double *h_effective_in,
     double *h_effective_reset,
     double beta_true, int len_gammas, double *gammas,
     int SAW_length_min, int SAW_length_max, int nmoves, int n_adaptations,
     int nruns_per_param_update,
     acq_maximizer_fn_t acq_maximizer_fn,
     gp_add_data_fn_t add_data_to_gp_fn,
     gp_maximizer_fn_t maximize_gp_fn,
     gp_posterior_mean_fn_t gp_posterior_mean_fn,
     gp_maximizer_local_fn_t gp_maximizer_local_fn,
     sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn,
     int strategy,
     double *E_samples, double *E_proposeds,
     double *log_f_fwds, double *log_f_revs,
     double *MH_ratios, int *SAW_lengths, double *used_gammas, double *rewards,
     double *policy, int policy_size)
{
    IMInputModel in_model;
    IMSimParams sim_params;
    IMMCMCDiagnostics out_data;

    build_input_model_ctypes
        (num_nodes, num_edges, J, h,
         nbrs_lists, nbrs_list_lens,
         incident_edge_lists, incident_edge_list_lens, edges,
         &in_model);
    populate_simulation_parameters_ctypes
        (h, S_ref, S_in, S_reset, h_effective_in, h_effective_reset,
         beta_true, len_gammas, gammas,
         SAW_length_min, SAW_length_max,
         nmoves, n_adaptations, nruns_per_param_update, num_nodes,
         &sim_params);

    int nactions = (SAW_length_max - SAW_length_min + 1) * len_gammas;

    set_up_common_IMRun_output_args_ctypes
        (E_samples, E_proposeds, log_f_fwds, log_f_revs,
         MH_ratios, SAW_lengths, used_gammas, rewards,
         nactions, S_end,
         &out_data);

    sampleIM_ego((const IMInputModel*)&in_model,
                 (const IMSimParams*)&sim_params,
                  &out_data,
                  acq_maximizer_fn, add_data_to_gp_fn, maximize_gp_fn,
                  gp_posterior_mean_fn, gp_maximizer_local_fn, sample_gp_posterior_mean_fn, (size_t)strategy, policy, (size_t)policy_size);

    free_input_model(&in_model);
}

void build_input_model_ctypes
    (int num_nodes, int num_edges,
     double *J, double *h,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens,
     int *edges,
     IMInputModel *in_model)
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
    (double *h, int *S_ref, int *S_in, int *S_reset,
     double *h_effective_in, double *h_effective_reset,
     double beta_true, int len_gammas, double *gammas,
     int SAW_length_min, int SAW_length_max,
     int nmoves, int n_adaptations, int nruns_per_param_update, int num_nodes,
     IMSimParams *sim_params)
{
    sim_params->S_ref = (const int*)S_ref;
    sim_params->S_in = (const int*)S_in;
    sim_params->S_reset = (const int*)S_reset;

    sim_params->h_effective_in = (const double*)h_effective_in;
    sim_params->h_effective_reset = (const double*)h_effective_reset;

    sim_params->E_in =
        compute_E_in(sim_params->S_in, sim_params->h_effective_in,
                     h, (int)num_nodes);

    sim_params->beta_true = beta_true;

    sim_params->gammas = (const double*)gammas;
    sim_params->n_gammas = (int)len_gammas;

    sim_params->SAW_length_min = SAW_length_min;
    sim_params->SAW_length_max = SAW_length_max;
    sim_params->n_moves = (int)nmoves;

    sim_params->n_adaptations = n_adaptations;

    sim_params->nruns_per_param_update = (int)nruns_per_param_update;
}

void set_up_common_IMRun_output_args_ctypes
    (double *E_samples, double *E_proposeds,
     double *log_f_fwds, double *log_f_revs,
     double *MH_ratios, int *SAW_lengths, double *used_gammas,
     double *rewards, int nactions, int *S_end,
     IMMCMCDiagnostics *out_data)
{
    out_data->E_samples = E_samples;
    out_data->E_proposed = E_proposeds;
    out_data->log_f_fwd = log_f_fwds;
    out_data->log_f_rev = log_f_revs;
    out_data->MH_ratio = MH_ratios;

    out_data->SAW_length = (int *)SAW_lengths;
    out_data->gammas = (double *)used_gammas;

    out_data->rewards = rewards;

    out_data->n_actions = (int)nactions;
    out_data->S_end = S_end;
}

