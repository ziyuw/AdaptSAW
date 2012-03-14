#ifndef __SARD_EGO_H_
#define __SARD_EGO_H_

#include "SARD.h"

#include "ego_callbacks.h"

void sampleSARD_ego(const SARDInputModel* p_in_model,
                    SARDSimParams* p_sim_params,
                    SARDMCMCDiagnostics* p_out_data,
                    acq_maximizer_fn_t acq_maximizer_fn,
                    gp_add_data_fn_t add_data_to_gp_fn,
                    gp_maximizer_fn_t maximize_gp_fn,
                    gp_posterior_mean_fn_t gp_posterior_mean_fn,
		    gp_maximizer_local_fn_t gp_maximizer_local_fn,
		    sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn,
		    return_sample_fn_t return_sample_fn, size_t return_samples,
                    size_t strategy, double *policy, size_t policy_size);

#endif // __SARD_EGO_H_

