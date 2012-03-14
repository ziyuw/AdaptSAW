#include "IM_ego.h"
#include "IMUtils.h"
#include <math.h>
#include <stdio.h>
#include "objectives.h"
#include "sampling.h"
#include <float.h>

#include <mt19937ar.h>

void 
sampleIM_ego(const IMInputModel * p_in_model, const IMSimParams * p_sim_params,
             IMMCMCDiagnostics * p_out_data,
             acq_maximizer_fn_t acq_maximizer_fn,
             gp_add_data_fn_t add_data_to_gp_fn,
             gp_maximizer_fn_t maximize_gp_fn,
             gp_posterior_mean_fn_t gp_posterior_mean_fn,
	     gp_maximizer_local_fn_t gp_maximizer_local_fn,
	     sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn,
	     size_t strategy,
             double *policy, size_t policy_size)
{
    
    double          E_initial;
    SAWData         saw_data;
    IMState         im_state;

    initIMState(&im_state, p_in_model, p_sim_params);

    // Set up the heaps, etc.
    initSAWData(&saw_data, p_in_model, p_sim_params);

    int nruns_until_next_param_update = p_sim_params->nruns_per_param_update;
    int first_run_with_curr_param_idx = 0;

    // Set up the bounds for EGO. The bounds are made up of ndim*2 numbers:
    // [bounds[2i], bounds[2i+1]] are the [lb, ub] of the ith dimension.
    // There are two parameters, SAW length and gamma, so there are two
    // dimensions. Each dimension has two bounds. The total size is 4.
    // The first dimension is SAW length. The second dimension is gamma.
    const size_t n_amcmc_dims = 2;
    double *bounds = malloc(sizeof(double) * n_amcmc_dims * 2);
    bounds[0] = p_sim_params->SAW_length_min;
    bounds[1] = p_sim_params->SAW_length_max;
    bounds[2] = p_sim_params->gammas[0];
    bounds[3] = p_sim_params->gammas[p_sim_params->n_gammas-1];

    double *amcmc_params = malloc(sizeof(double) * n_amcmc_dims);

    // Select initial parameters - this must be done without the help of EGO
    // because the Gaussian process contains no information at all. The initial
    // parameters can either be picked randomly or set to something reasonable.
    // acq_maximizer_fn(bounds, amcmc_params, n_amcmc_dims);
    amcmc_params[0] = (bounds[1] + bounds[0]) / 2;
    amcmc_params[1] = (bounds[3] + bounds[2]) / 2;

    saw_data.SAW_length = (int)round(amcmc_params[0]);
    saw_data.gamma_down = amcmc_params[1];
    saw_data.gamma_up = amcmc_params[1];

    objective obj_fn = ACF;

    double target_avg_mh_ratio = 0.4;
    double scale_mult = 1 / max(1 - target_avg_mh_ratio, target_avg_mh_ratio);

    // the first adaptation is already done (and is hardcoded, since the GP has
    // no information at the start of the adaptation stage).
    int n_adaptations = p_sim_params->n_adaptations - 1;
    int n_adaptations_done = 0;
    int adaptation_finished = 0;

    // size_t policy_size = 1000;
    // double *policy = malloc(sizeof(double) * n_amcmc_dims * policy_size);

    for(int curr_move = 0; curr_move < p_sim_params->n_moves; curr_move++)
    {
	if (curr_move%1000 == 0) {
	    printf("Iteration %d Finished.\n", curr_move);
	}
	
        if (!adaptation_finished &&
            nruns_until_next_param_update == 0 && curr_move > 0)
        {
            double reward = -1.0;

            if (obj_fn == ACF)
            {
                reward =
                    compute_robust_avg_abs_acf_est
                        (p_out_data->E_samples,
                         p_sim_params->nruns_per_param_update,
                         first_run_with_curr_param_idx);
            }
            else if (obj_fn == MH_RATIO)
            {
                double *mh_ratio_window_start =
                    p_out_data->MH_ratio+first_run_with_curr_param_idx;
                int mh_ratio_window_size =
                    p_sim_params->nruns_per_param_update;

                double avg_mh_ratio =
                    compute_mean(mh_ratio_window_start,
                            mh_ratio_window_size);

                reward = 1.0 -
                    scale_mult *
                    fabs(target_avg_mh_ratio - avg_mh_ratio);
            }
            else if (obj_fn == MH_RATIO_AVG_ABS_DIFF)
            {
                double *mh_ratio_window_start =
                    p_out_data->MH_ratio+first_run_with_curr_param_idx;
                int mh_ratio_window_size =
                    p_sim_params->nruns_per_param_update;

                double mh_ratio_avg_abs_diff = 0.0;

                for (int i = 0; i < mh_ratio_window_size; ++i)
                {
                    mh_ratio_avg_abs_diff +=
                        fabs(mh_ratio_window_start[i] -
                                target_avg_mh_ratio);
                }

                mh_ratio_avg_abs_diff /= mh_ratio_window_size;

                reward = 1.0 - scale_mult * mh_ratio_avg_abs_diff;
            }
            else
            {
                printf("Invalid objective function specified\n");
            }

            // Update the Gaussian process
            add_data_to_gp_fn(amcmc_params, n_amcmc_dims, reward);

            // IMPL: Check whether the best possible parameter settings have
            // improved much and if they haven't, then stop adapting. I'm not
            // sure exactly how to do this..

            for (int i = 0; i < p_sim_params->nruns_per_param_update; ++i)
            {
                p_out_data->rewards[i+first_run_with_curr_param_idx] =
                    reward;
            }
                
            first_run_with_curr_param_idx += p_sim_params->nruns_per_param_update;

            if (n_adaptations > 0)
            {
                // Select next set of parameters
                acq_maximizer_fn(bounds, amcmc_params, n_amcmc_dims);
                saw_data.SAW_length = (int)round(amcmc_params[0]);
                saw_data.gamma_down = amcmc_params[1];
                saw_data.gamma_up = amcmc_params[1];

                --n_adaptations;
                ++n_adaptations_done;
            }
            else if (n_adaptations == 0)
            {
                // adaptation is done - construct probability distribution over
                // parameters and draw parameters from it
                policy_size = build_sir_policy(policy, policy_size, bounds, strategy,
                                 gp_posterior_mean_fn, acq_maximizer_fn, gp_maximizer_local_fn, sample_gp_posterior_mean_fn);

                size_t param_idx = sample_uniform(policy_size);
                saw_data.SAW_length = (int)round(policy[2*param_idx]);
                saw_data.gamma_down = policy[2*param_idx+1];
                saw_data.gamma_up = policy[2*param_idx+1];

                // adaptation is done - reset IM state
                setIMState(&im_state, p_in_model,
                           (const int*)p_sim_params->S_reset,
                           (const double*)p_sim_params->h_effective_reset);

                adaptation_finished = 1;
            }

            nruns_until_next_param_update = p_sim_params->nruns_per_param_update;
        }

        if (!adaptation_finished)
        {
            --nruns_until_next_param_update;
        }
        else
        {
            size_t param_idx = sample_uniform(policy_size);
            saw_data.SAW_length = (int)round(policy[2*param_idx]);
            saw_data.gamma_down = policy[2*param_idx+1];
            saw_data.gamma_up = policy[2*param_idx+1];
        }

        E_initial = im_state.E;
        proposeIMMove(&saw_data, &im_state, p_in_model, p_sim_params);

        // Compute MH_ratio
        double MH_ratio = p_sim_params->beta_true * (E_initial - im_state.E);
        MH_ratio += (saw_data.log_f_rev - saw_data.log_f_fwd);
        MH_ratio = exp(MH_ratio);
        MH_ratio = MH_ratio > 1.0 ? 1.0 : MH_ratio;

        storeDiagnosticsData(p_out_data, curr_move, E_initial, im_state.E,
                             saw_data.log_f_fwd, saw_data.log_f_rev, MH_ratio,
                             saw_data.SAW_length, saw_data.gamma_up);

        // If rejection, restore state:
        if (genrand_real2() > MH_ratio)
        {
            undoIMMove(&saw_data, &im_state, p_in_model, p_sim_params);
        }
    }

    storeEndState(p_out_data, im_state.S, p_in_model->num_nodes);

    freeIMState(&im_state, p_in_model, p_sim_params);
    freeSAWData(&saw_data, p_in_model, p_sim_params);

    // free(policy);
    free(amcmc_params);
    free(bounds);
}

size_t build_sir_policy(double *policy, size_t policy_size,
                      double *bounds, size_t strategy,
                      gp_posterior_mean_fn_t gp_posterior_mean_fn,
		      acq_maximizer_fn_t acq_maximizer_fn,
		      gp_maximizer_local_fn_t gp_maximizer_local_fn,
		      sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn)
{
    printf("build_sir_policy(): START with strategy: %d\n", strategy);
    
    // strategy: 	0 Nim's Resampling
    //			1 l-bfgs-b + resampling
    //			2 threshholding
    //			3 sampling + direct
    //			4 sampling + multiple l-bfgs-b
    
    if (strategy == 0) {
	double SAW_length_min = bounds[0];
	double SAW_length_max = bounds[1];
	double gamma_min = bounds[2];
	double gamma_max = bounds[3];

	size_t n_gammas = 100;

	// Draw sample of points
	double *grid_points = \
	    make_deterministic_grid_points
		(SAW_length_max, SAW_length_min, gamma_max, gamma_min, n_gammas);

	size_t n_sawlengths = (int)(SAW_length_max - SAW_length_min + 1);
	size_t n_points = n_sawlengths * n_gammas;

	// Compute w(x) = f(\mu(x)) for each point
	double *w = malloc(sizeof(double) * n_points);

	for (size_t i = 0; i < n_points; ++i)
	{
	    w[i] = gp_posterior_mean_fn(grid_points+(2*i), 2);
	    w[i] = exp(10*w[i]);
	}

	// Compute normalization constant
	double Z = 0.0;
	for (size_t i = 0; i < n_points; ++i)
	{
	    Z += w[i];
	}

	// Normalize weights
	for (size_t i = 0; i < n_points; ++i)
	{
	    w[i] /= Z;
	}

	// Resample points according to weights to construct policy
	size_t *policy_idx = malloc(sizeof(size_t) * policy_size);
	kitagawa_resampling(w, n_points, policy_idx, policy_size);
	for (size_t i = 0; i < policy_size; ++i)
	{
	    policy[2*i] = grid_points[2*policy_idx[i]];
	    policy[2*i+1] = grid_points[2*policy_idx[i]+1];
	}

	free(policy_idx);
	free(w);
	free(grid_points);

    }else if (strategy >= 1 ) {
	// Options: 	0 l-bfgs-b + resampling
	//		1 threshholding
	//		2 sampling + direct
	//		3 sampling + multiple l-bfgs-b
	
	size_t options = strategy - 1;
	
	const size_t n_amcmc_dims = 2;
	double *amcmc_params, *amcmc_params_all, *w;
	size_t counter = 0;
	amcmc_params = malloc(sizeof(double) * (n_amcmc_dims+1));
	amcmc_params_all = malloc(sizeof(double) * n_amcmc_dims * policy_size);
	w = malloc(sizeof(double) * policy_size);
	size_t length = 0;
	
	// set w[i] to 0 for all i
	for (size_t i = 0 ; i < policy_size; i++) {
	    w[i] = 0.0;
	}
	
	size_t num_loop = policy_size;
	double max_val = -100.0;
	double * amcmc_params_max;
	
	if (options == 2) {
	    num_loop = 1;
	}
	
	if (options == 3) {
	    num_loop = 100;
	    amcmc_params_max = malloc(sizeof(double) * (n_amcmc_dims+1));
	}
	
	while (counter < num_loop) {
	    gp_maximizer_local_fn(bounds, amcmc_params, n_amcmc_dims);
	    
	    if (options == 2) {
		sample_gp_posterior_mean_fn(amcmc_params, n_amcmc_dims, bounds, policy, policy_size);
	    } else if (options < 2) {
		check_redundancy(amcmc_params_all, amcmc_params, w, &length); 
		printf("length is now: %d\n", length);
	    } else if (options == 3) {
		if (amcmc_params[2] > max_val) {
		    max_val = amcmc_params[2];
		    for (size_t k = 0 ; k <= n_amcmc_dims ; k++) {
			amcmc_params_max[k] = amcmc_params[k];
		    }
		}
	    }
	    
	    counter = counter + 1;
	}
	
	if (options == 3) {
	    sample_gp_posterior_mean_fn(amcmc_params_max, n_amcmc_dims, bounds, policy, policy_size);
	    free(amcmc_params_max);
	}
	
	if (options == 1) {
	    thresholding(amcmc_params_all, w, &length);
	    policy_size = length;
	    for (size_t i = 0; i < policy_size; ++i)
	    {
		policy[2*i] = amcmc_params_all[2*i];
		policy[2*i+1] = amcmc_params_all[2*i+1];
	    }
	} else if (options == 0) {
	    // Compute normalization constant
	    double Z = 0.0;
	    for (size_t i = 0; i < policy_size; ++i)
	    {
		Z += w[i];
	    }

	    // Normalize weights
	    for (size_t i = 0; i < policy_size; ++i)
	    {
		w[i] /= Z;
	    }
	    
	    size_t *policy_idx = malloc(sizeof(size_t) * policy_size);
	    kitagawa_resampling(w, length, policy_idx, policy_size);
	    
	    for (size_t i = 0; i < policy_size; ++i)
	    {
		policy[2*i] = amcmc_params_all[2*policy_idx[i]];
		policy[2*i+1] = amcmc_params_all[2*policy_idx[i]+1];
	    }
	    free(policy_idx);
	}
	
	free(w);
	free(amcmc_params_all);
	free(amcmc_params);
    }
    
    printf("build_sir_policy(): FINISH\n");
    return policy_size;
}

void check_redundancy(double* existing_pts, double *new_pt, double *w, size_t* length) {
    // First round new_pt
    new_pt[0] = floor(new_pt[0]);
    new_pt[1] = floor(new_pt[1]*100.0)/100.0;
    
    for (size_t i = 0 ; i < *length ; i++) {
	if (new_pt[0] == existing_pts[i*2] && new_pt[1] == existing_pts[i*2+1]) {
	    return;
	}
    }
    
    existing_pts[*length*2] = new_pt[0];
    existing_pts[*length*2+1] = new_pt[1];
    w[*length] = new_pt[2];
    *length = *length + 1;
}

void thresholding(double* all_pts, double *w, size_t* length) {
    double max_val = -10.0;
    for (size_t i = 0 ; i < *length ; i++) {
	if (w[i] > max_val) {
	    max_val = w[i];
	}
    }
    size_t counter = 0;
    
    for (size_t i = 0 ; i < *length ; i++) {
	if (w[i] > 0.8*max_val) {
	    all_pts[2*counter] = all_pts[i*2];
	    all_pts[2*counter+1] = all_pts[i*2+1];
	    w[counter] = w[i];
	    counter = counter + 1;
	}
    }
    
    for (size_t i = counter ; i < *length ; i++) {
	all_pts[2*counter] = 0;
	all_pts[2*counter+1] = 0;
	w[i] = 0.0;
    }
    
    *length = counter;
}

double* make_deterministic_grid_points
    (double SAW_length_max, double SAW_length_min,
     double gamma_max, double gamma_min, size_t n_gammas)
{
    size_t n_amcmc_dims = 2;

    size_t n_sawlengths = (int)(SAW_length_max - SAW_length_min + 1);
    size_t n_points = n_sawlengths * n_gammas;
    double *points = malloc(sizeof(double) * n_amcmc_dims * n_points);

    double gamma_incr = (gamma_max - gamma_min) / (double)n_gammas;

    for (size_t i = 0; i < n_gammas; ++i)
    {
        for (size_t j = 0; j < n_sawlengths; ++j)
        {
            size_t param_idx = i*n_sawlengths + j;
            points[2*param_idx] = SAW_length_min + (double)j;
            points[2*param_idx+1] = gamma_min + (i * gamma_incr);
        }
    }

    return points;
}

