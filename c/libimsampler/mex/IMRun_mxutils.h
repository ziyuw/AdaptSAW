#ifndef __IMRUN_MXUTILS_H_
#define __IMRUN_MXUTILS_H_

#include "mex.h"
#include "IM.h"

void set_up_common_IMRun_output_args_mex
    (int nlhs, mxArray *plhs[], IMMCMCDiagnostics *out_data,
     const int n_moves, const int n_actions);

void populate_simulation_parameters_mex
    (int nrhs, const mxArray *prhs[], IMSimParams *sim_params, int num_nodes);
void build_input_model_mex
    (int nrhs, const mxArray *prhs[], IMInputModel *in_model);

#endif // __IMRUN_MXUTILS_H_

