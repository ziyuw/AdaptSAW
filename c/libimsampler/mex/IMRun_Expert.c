/* Matlab API to the IM function */

/* 3.Jan.10: DO NOT call this function on "dense" graphs as it uses
   heaps! */

#include "mex.h"
#include "IMOutData.h"
#include "IM.h"
#include "IMRun_mxutils.h"
#include "IM_expert.h"
#include "IMUtils.h"

// CURRENT compilation command: Use the Makefile target IMRun_Expert.mex

void mexFunction(  
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    /* USAGE:
       [E_samples_vector, E_proposed_vector,
        log_f_fwd_vector, log_f_rev_vector,
        MH_ratio_vector, SAW_length_vector,
        gains_vector, rewards_vector] =
       IMRun_Expert(nbrs_list_i, incident_edge_list_i, edges_i, J, h,
                    int32(S_ref), int32(S_in), h_effective_in, beta_true,
                    gammas, int32(SAW_length_min), int32(SAW_length_max),
                    int32(n_moves), int32(nruns_per_param_update));
       Rewards vector is just a vector of undefined values.
     */

    IMInputModel in_model;
    IMSimParams sim_params;
    IMMCMCDiagnostics out_data;

    build_input_model_mex(nrhs, prhs, &in_model);
    populate_simulation_parameters_mex(nrhs, prhs, &sim_params, in_model.num_nodes);

    int n_SAW_lengths =
        sim_params.SAW_length_max - sim_params.SAW_length_min + 1;
    int n_actions = n_SAW_lengths * sim_params.n_gammas;

    set_up_common_IMRun_output_args_mex
        (nlhs, plhs, &out_data, sim_params.n_moves, n_actions);

    sampleIM_expert((const IMInputModel*)&in_model,
                    (const IMSimParams*)&sim_params,
                    &out_data);

    free_input_model(&in_model);
}

