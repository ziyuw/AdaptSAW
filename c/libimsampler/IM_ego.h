#ifndef __IM_EGO_H_
#define __IM_EGO_H_

#include "ego_callbacks.h"
#include "IM.h"

void sampleIM_ego(const IMInputModel* p_in_model,
                  const IMSimParams* p_sim_params,
                  IMMCMCDiagnostics * p_out_data,
                  acq_maximizer_fn_t acq_maximizer_fn,
                  gp_add_data_fn_t add_data_to_gp_fn,
                  gp_maximizer_fn_t maximize_gp_fn,
                  gp_posterior_mean_fn_t gp_posterior_mean_fn,
		  gp_maximizer_local_fn_t gp_maximizer_local_fn,
		  sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn,
		  size_t strategy,
                  double *policy, size_t policy_size);

size_t build_sir_policy(double *policy, size_t policy_size,
                      double *bounds, size_t strategy,
                      gp_posterior_mean_fn_t gp_posterior_mean_fn,
		      acq_maximizer_fn_t acq_maximizer_fn,
		      gp_maximizer_local_fn_t gp_maximizer_local_fn,
		      sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn);

double*
make_deterministic_grid_points
    (double SAW_length_max, double SAW_length_min,
     double gamma_max, double gamma_min, size_t n_gammas);

#endif // __IM_EGO_H_

