#include "IMRun_mxutils.h"
#include "GraphNodes.h"
#include "GraphNodes_mex.h"
#include "IMUtils.h"

// Common output parameters are:
// [E_samples_vector, E_proposed_vector,
//  log_f_fwd_vector, log_f_rev_vector,
//  MH_ratio_vector, SAW_length_vector,
//  log_gains_vector, rewards_vector]
void set_up_common_IMRun_output_args_mex
    (int nlhs, mxArray *plhs[], IMMCMCDiagnostics *out_data,
     const int n_moves, const int n_actions)
{
    /* Allocate the output data and set the pointers: */
    plhs[0] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->E_samples = mxGetPr( plhs[0] );

    plhs[1] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->E_proposed = mxGetPr( plhs[1] );

    plhs[2] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->log_f_fwd = mxGetPr( plhs[2] );

    plhs[3] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->log_f_rev = mxGetPr( plhs[3] );

    plhs[4] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->MH_ratio = mxGetPr( plhs[4] );

    {
        int *dims = (int*)malloc(2*sizeof(int));
        dims[0] = n_moves;
        dims[1] = 1;

        plhs[5] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);

        free(dims);
    }
    out_data->SAW_length = (int *)mxGetPr( plhs[5] );

    plhs[6] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->gammas = mxGetPr( plhs[6] );

    out_data->n_actions = n_actions;

    plhs[7] = mxCreateDoubleMatrix( n_moves, 1, mxREAL );
    out_data->rewards = mxGetPr( plhs[7] );
}

// Common input parameters are:
// IMRun_Exp3(nbrs_list_i, incident_edge_list_i, edges_i, J, h,
//            int32(S_ref), int32(S_in), h_effective_in, beta_true,
//            gammas, int32(SAW_length_min), int32(SAW_length_max),
//            int32(n_moves), int32(nruns_per_param_update));
void populate_simulation_parameters_mex
    (int nrhs, const mxArray *prhs[], IMSimParams *sim_params, int num_nodes)
{
    const double *h = mxGetPr(prhs[4]);
    sim_params->S_ref = (const int*)mxGetData(prhs[5]);
    sim_params->S_in = (const int*)mxGetData(prhs[6]);
    sim_params->h_effective_in = (const double*)mxGetPr(prhs[7]);

    sim_params->E_in =
        compute_E_in(sim_params->S_in, sim_params->h_effective_in,
                     h, num_nodes);

    /* Note my alternative to mxGetScalar... */
    sim_params->beta_true = *( mxGetPr(prhs[8]) );

    sim_params->gammas = (const double*)mxGetPr(prhs[9]);
    sim_params->n_gammas = (int)mxGetM(prhs[9]);

    sim_params->SAW_length_min = *( (const int*)mxGetData(prhs[10]) );
    sim_params->SAW_length_max = *( (const int*)mxGetData(prhs[11]) );
    sim_params->n_moves = *( (const int*)mxGetData(prhs[12]) );

    sim_params->n_adaptations = *( (const int*)mxGetData(prhs[13]) );
    sim_params->nruns_per_param_update = *( (const int*)mxGetData(prhs[14]) );
}

void build_input_model_mex
    (int nrhs, const mxArray *prhs[], IMInputModel *in_model)
{
    /* Pointer to a 2D matrix, stored as a 1D array Fortran style */
    const int* p_edges = (const int*)mxGetData(prhs[2]);

    /* Pointers to real vectors */
    const double *J = mxGetPr(prhs[3]);
    const double *h = mxGetPr(prhs[4]);

    int num_nodes = mxGetM(prhs[5]);
    int num_edges = mxGetM(prhs[2]);

    /* ptrs to mxArrays inputted; these 2 are cell arrays */
    const mxArray *mx_p_nbrs_list = prhs[0];
    const mxArray *mx_p_incident_edge_list = prhs[1];

    GraphNodes* graph_nodes =
        build_graphnodes_mex(mx_p_nbrs_list, mx_p_incident_edge_list, num_nodes);
    int** graph_edges = build_graphedges_mex(p_edges, num_edges);

    in_model->num_nodes = num_nodes;
    in_model->num_edges = num_edges;
    in_model->graph_nodes = (const GraphNodes*)graph_nodes;
    in_model->graph_edges = (const int**)graph_edges;
    in_model->J = J;
    in_model->h = h;
}

