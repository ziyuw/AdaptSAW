#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<stdbool.h>
#include<mt19937ar.h>
#include "GraphNodes.h"
#include "SARDOutData.h"
#include "SAWHeaps.h"
#include "SARD.h"
#include "SARDUtils.h"

#include "logspace.h"
#include "utils.h"
#include "sampling.h"

#include "ego_callbacks.h"
#include "objectives.h"

/* 
   NOTE: Unconstrained sampling on the "full" state space given an
   input state. For ergodicity need to be able to do single bit flips.

   3.Dec.10: This is very similar in spirit to V4 but a few differences:

   1) Separate gamma_hi and gamma_lo are input

   2) Performs iterated SAWs instead of single-flips

   3) Has the option of calculating the overlap of 2 independent realization (important to check mixing!)
   -- Yea I took this out. - Nim

   
   TODO: 
   
*/

void assign_sard_params
    (SARDState* p_sard_state, SAWData *p_saw_data, const SARDInputModel *p_in_model,
     SARDSimParams *p_sim_params, double *sard_params, int normalize_mix_p);
int check_violated_sard_constraints(double *amcmc_params, int n_amcmc_dims);

void build_sir_policy(double *policy, size_t policy_size,
                      double *bounds, size_t ndims,
		      size_t strategy,
                      gp_posterior_mean_fn_t gp_posterior_mean_fn,
		      acq_maximizer_fn_t acq_maximizer_fn,
		      gp_maximizer_local_fn_t gp_maximizer_local_fn,
		      sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn);

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
                    size_t strategy, double *policy, size_t policy_size)
{
    SAWData* p_saw_data = (SAWData*)malloc(sizeof(SAWData));
    SARDState* p_sard_state = (SARDState*)malloc(sizeof(SARDState));

    initSARDState(p_sard_state, p_in_model, p_sim_params);
    // Set up the heaps, etc. 
    initSAWData(p_saw_data, p_in_model, p_sim_params);

    int nruns_until_next_param_update = p_sim_params->nruns_per_param_update;
    int first_run_with_curr_param_idx = 0;

    double *bounds, *amcmc_params;
    const size_t n_amcmc_dims = 8;
    {
        // Set up the bounds for EGO. The bounds are made up of ndim*2 numbers:
        // [bounds[2i], bounds[2i+1]] are the [lb, ub] of the ith dimension.
        // There are eight parameters, iters_per_move, SAW_length_[min, max],
        // P_LL, P_LH, P_HL, gamma_H and gamma_L, so there are eight
        // dimensions. Each dimension has two bounds. The total size is 16.
        bounds = malloc(sizeof(double) * n_amcmc_dims * 2);
        bounds[0] = (double)p_sim_params->iters_per_move_min;
        bounds[1] = (double)p_sim_params->iters_per_move_max;
        bounds[2] = (double)p_sim_params->min_SAW_length_bound_min;
        bounds[3] = (double)p_sim_params->min_SAW_length_bound_max;
        bounds[4] = (double)p_sim_params->max_SAW_length_bound_min;
        bounds[5] = (double)p_sim_params->max_SAW_length_bound_max;
        bounds[6] = 0.0; // P_LL
        bounds[7] = 1.0; // P_LL
        bounds[8] = 0.0; // P_LH
        bounds[9] = 1.0; // P_LH
        bounds[10] = 0.0; // P_HL
        bounds[11] = 1.0; // P_HL
        bounds[12] = p_sim_params->gamma_hi_min;
        bounds[13] = p_sim_params->gamma_hi_max;
        bounds[14] = p_sim_params->gamma_lo_min;
        bounds[15] = p_sim_params->gamma_lo_max;

        // Select initial parameters - this must be done without the help of EGO
        // because the Gaussian process contains no information at all. The initial
        // parameters can either be picked randomly or set to something reasonable.
        // acq_maximizer_fn(bounds, amcmc_params, n_amcmc_dims);
        amcmc_params = malloc(sizeof(double) * n_amcmc_dims);
        for (int i = 0; i < n_amcmc_dims; ++i)
        {
            int idx = i * 2;
            amcmc_params[i] = (bounds[idx + 1] + bounds[idx]) / 2;
        }

        assign_sard_params(p_sard_state, p_saw_data, p_in_model,
                           p_sim_params, amcmc_params, 1);

        int constraints_violated =
            check_violated_sard_constraints(amcmc_params, n_amcmc_dims);
        if (constraints_violated)
        {
            printf("Initial SARD_ego parameters violate constraints.\n");
        }
    }

    objective obj_fn = ACF;

    // the first adaptation is already done (and is hardcoded, since the GP has
    // no information at the start of the adaptation stage).
    int n_adaptations = p_sim_params->n_adaptations - 1;
    int n_adaptations_done = 0;
    int adaptation_finished = 0;

    for (int curr_move=0; curr_move<p_sim_params->n_moves; curr_move++)
    {
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
            else
            {
                printf("Invalid objective function specified\n");
            }

	    printf("Reward: %f\n", reward);
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
                int num_times_violated = 0;

                while (true)
                {
                    acq_maximizer_fn(bounds, amcmc_params, n_amcmc_dims);

                    int constraints_violated =
                        check_violated_sard_constraints(amcmc_params, n_amcmc_dims);
                    if (!constraints_violated)
                    {
                        break;
                    }

                    ++num_times_violated;
                    printf("Constraints violated %d times\n",
                            num_times_violated);

                    add_data_to_gp_fn(amcmc_params, n_amcmc_dims, 0.0);
                }

                assign_sard_params(p_sard_state, p_saw_data, p_in_model,
                                   p_sim_params, amcmc_params, 1);

                --n_adaptations;
                ++n_adaptations_done;
            }
            else if (n_adaptations == 0)
            {
                // adaptation is done - construct probability distribution over
                // parameters and draw parameters from it
                printf("START: Building SIR policy\n");
                build_sir_policy(policy, policy_size, bounds, n_amcmc_dims, strategy,
                                 gp_posterior_mean_fn, acq_maximizer_fn, gp_maximizer_local_fn,
				 sample_gp_posterior_mean_fn);
                printf("STOP: Building SIR policy\n");
                size_t param_idx = sample_uniform(policy_size);
                
                assign_sard_params
                    (p_sard_state, p_saw_data, p_in_model, p_sim_params,
                     policy+(n_amcmc_dims*param_idx), 1);  

                
		int reset = 0;
		if (reset) {
		    // adaptation is done - reset SARD state
		    printf("START: Resetting SARD state\n");
		    setSARDState(p_sard_state, p_in_model,
				(const int*)p_sim_params->S_reset,
				(const double*)p_sim_params->h_effective_reset);
		    printf("STOP: Resetting SARD state\n");
		}
                adaptation_finished = 1;
            }
            else
            {
                printf("Invalid value for n_adaptations\n");
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
            assign_sard_params(p_sard_state, p_saw_data, 
			       p_in_model, p_sim_params, 
			       policy + (n_amcmc_dims*param_idx), 1);
        }

        // Select move type:
        {
            double R=genrand_real2();
            if (R <= p_sim_params->P_LL)
            {
                p_saw_data->move_type=LL;
            }
            else if (R <= p_sim_params->P_LL+p_sim_params->P_LH)
            {
                p_saw_data->move_type=LH;
            }
            else
            {
                p_saw_data->move_type=HL;
            }
        }

        double E_initial= p_sard_state->E;
        proposeSAWMove(p_saw_data, p_sard_state, p_in_model, p_sim_params);

        double log_f_fwd, log_f_rev;
        evalLogFs(&log_f_fwd, &log_f_rev, p_saw_data, p_sim_params);

        double MH_ratio = p_sim_params->beta_true*(E_initial - p_sard_state->E);
        MH_ratio += (log_f_rev - log_f_fwd);
        MH_ratio = exp(MH_ratio);
        MH_ratio = MH_ratio >1.0 ? 1.0 : MH_ratio;

        double reward = 0.0;

	if (adaptation_finished && return_samples == 1) {
	    // return_samples specifies whether samples are to be returned
	    // 0: do not return
	    // 1: return samples after adaptation is finished
	    // 2: return all samples	    
	    return_sample_fn(p_sard_state->S, p_in_model->num_nodes);
	}
	
        storeDiagnosticsData(p_out_data, p_sim_params, curr_move,
                             E_initial, p_sard_state->E,
                             log_f_fwd, log_f_rev, MH_ratio,
                             (int)p_saw_data->move_type,
                             p_saw_data->sigma_lengths,
                             p_saw_data->rho_lengths,
                             p_sim_params->iters_per_move,
                             p_sim_params->gamma_hi, p_sim_params->gamma_lo,
                             p_sim_params->P_LL, p_sim_params->P_LH,
                             p_sim_params->P_HL, reward);

        // If rejection, restore state:
        if (genrand_real2() > MH_ratio)
        { 
            undoSAWMove(p_saw_data, p_sard_state, p_in_model, p_sim_params);
        }
        
        if ((curr_move+1)%1000 == 0) {
	    printf("Iteration %d Finished.\n", curr_move+1);
	}
    } // end moves loop

    printf("X\n");

    storeEndState(p_out_data, p_sard_state->S, p_in_model->num_nodes);

    printf("Y\n");

    freeSARDState(p_sard_state, p_in_model, p_sim_params);
    freeSAWData(p_saw_data, p_in_model, p_sim_params);  
    free(p_sard_state);
    free(p_saw_data);
    free(amcmc_params);
    free(bounds);
    printf("Z\n");
}

void assign_sard_params
    (SARDState* p_sard_state, SAWData *p_saw_data, const SARDInputModel *p_in_model,
     SARDSimParams *p_sim_params, double *sard_params, int normalize_mix_p)
{
    p_sim_params->iters_per_move = (int)round(sard_params[0]);
    p_sim_params->SAW_length_min = (int)round(sard_params[1]);
    p_sim_params->SAW_length_max = (int)round(sard_params[1]+sard_params[2]);
    
    p_sim_params->P_LL = sard_params[3];
    p_sim_params->P_LH = sard_params[4];
    p_sim_params->P_HL = sard_params[5];
    
    p_sim_params->gamma_hi = sard_params[7]+sard_params[6];
    p_sim_params->gamma_lo = sard_params[7];
    
    if (!normalize_mix_p)
    {
	resetSAWHeaps(p_saw_data, p_in_model, p_sim_params);
        return;
    }

    // Pick P_LL, P_LH and P_HL, satisfying the constraint that
    // P_LL + P_LH + P_HL = 1 by normalizing.
    double mix_Z =
        exp(log_sum_exp(log_sum_exp(log(p_sim_params->P_LL),
                                    log(p_sim_params->P_HL)),
                        log(p_sim_params->P_LH)));

    p_sim_params->P_LL /= mix_Z;
    p_sim_params->P_LH /= mix_Z;
    p_sim_params->P_HL /= mix_Z;
    
    int debug = 0;
    if (debug) {
	printf("iters_per_move: %d\n", p_sim_params->iters_per_move);
	printf("length_min: %d\n", p_sim_params->SAW_length_min);
	printf("length_max: %d\n", p_sim_params->SAW_length_max);
	
	printf("p_sim_params->P_LL: %f\n", p_sim_params->P_LL);
	printf("p_sim_params->P_LH: %f\n", p_sim_params->P_LH);
	printf("p_sim_params->P_HL: %f\n", p_sim_params->P_HL);
	
	printf("p_sim_params->gamma_hi: %f\n", p_sim_params->gamma_hi);
	printf("p_sim_params->gamma_lo: %f\n", p_sim_params->gamma_lo);
    }
    
    resetSAWHeapsCurrent(p_sard_state, p_saw_data, p_in_model, p_sim_params);
}

int check_violated_sard_constraints(double *amcmc_params, int n_amcmc_dims)
{
// NOTE: This method is unnecessary anymore but I kept it to minimize the change necessary.

//     int SAW_length_min = (int)round(amcmc_params[1]);
//     int SAW_length_max = (int)round(amcmc_params[2]);
// 
//     if (SAW_length_min > SAW_length_max)
//     {
//         return 1;
//     }

    return 0;
}

void build_sir_policy(double *policy, size_t policy_size,
                      double *bounds, size_t ndims,
		      size_t strategy,
                      gp_posterior_mean_fn_t gp_posterior_mean_fn,
		      acq_maximizer_fn_t acq_maximizer_fn,
		      gp_maximizer_local_fn_t gp_maximizer_local_fn,
		      sample_gp_posterior_mean_fn_t sample_gp_posterior_mean_fn)
{
    // strategy: 	0 Nim's Resampling
    //			1 l-bfgs-b + M_resampling
    //			2 threshholding
    //			3 sampling + direct
    //			4 sampling + multiple l-bfgs-b
    
    printf("build_sir_policy(): START\n");
    
    if (strategy > 0) {
	const size_t n_amcmc_dims = 8;
	double *amcmc_params;
	double *amcmc_params_max;
	size_t counter = 0;
	amcmc_params = malloc(sizeof(double) * (n_amcmc_dims+1));
	
	size_t num_loops = policy_size;
	double max_val = -100.0;
	
	if (strategy == 1) {
	    num_loops = 100;
	    amcmc_params_max = malloc(sizeof(double) * (n_amcmc_dims));
	} else if (strategy == 2) {
	    num_loops = 2;
	    amcmc_params_max = malloc(sizeof(double) * (n_amcmc_dims));
	}
	
	while (counter < num_loops) {
	    gp_maximizer_local_fn(bounds, amcmc_params, n_amcmc_dims);
	    
	    if (strategy == 1 || strategy == 2) {
		if (amcmc_params[n_amcmc_dims] > max_val) {
		    max_val = amcmc_params[n_amcmc_dims];
		    for (size_t k = 0 ; k < n_amcmc_dims ; k++) {
			amcmc_params_max[k] = amcmc_params[k];
		    }
		}
		counter = counter + 1;
	    } else if (strategy == 3) {
		int constraints_violated = check_violated_sard_constraints(amcmc_params, n_amcmc_dims);
		if (!constraints_violated)
		{
		    amcmc_params[1] = round(amcmc_params[1]);
		    amcmc_params[2] = round(amcmc_params[2]);
		    
		    for (size_t d = 0; d < ndims; ++d)
		    {
			policy[counter*ndims+d] = amcmc_params[d];
		    }
		    counter = counter + 1;
		} else {
		    printf("Constraint Violated!\n");
		}    
	    }
	}
	
	if (strategy == 1 || strategy == 2) {
	    sample_gp_posterior_mean_fn(amcmc_params_max, n_amcmc_dims, bounds, policy, policy_size);
	    free(amcmc_params_max);
	}
    } else {
    
	// Do not put any sets of parameters that violate constraints into the
	// policy. It just makes things more complicated.

	// Compute number of points per dimension. This will probably result in way
	// too many points, so it would probably be a good idea to construct the
	// points in a more intelligent way.

	printf("1\n");
	size_t *npoints_per_dim = malloc(sizeof(size_t) * ndims);
	{
	    for (size_t i = 0; i < 3; ++i)
	    {
		npoints_per_dim[i] = (size_t)(bounds[2*i+1] - bounds[2*i] + 1);
	    }

	    for (size_t i = 3; i < ndims; ++i)
	    {
		npoints_per_dim[i] = 6;
	    }
	}

	printf("2\n");
	// Make grid of points deterministically.
	size_t n_points = 30000;
	double *grid_points =
	    make_uniform_grid_points_nd(bounds, npoints_per_dim, ndims, n_points);

	printf("3\n");
	// Compute w(x) = f(\mu(x)) for each point. Assign w(x) = 0 to the points
	// that violate constraints.
	double *w = malloc(sizeof(double) * n_points);

	printf("4\n");

	printf("computing gp posterior for %lu points\n", n_points);
	size_t n_violated_points = 0, n_valid_points = 0;
	for (size_t i = 0; i < n_points; ++i)
	{
	    double *curr_point = grid_points+(ndims*i);

	    if (check_violated_sard_constraints(curr_point, ndims))
	    {
		w[i] = 0.0;
		++n_violated_points;
		continue;
	    }

	    ++n_valid_points;

	    w[i] = gp_posterior_mean_fn(curr_point, ndims);
	    w[i] = exp(10*w[i]);

	    if (i % 1000 == 0)
	    {
		printf("n_violated = %lu, n_valid = %lu, total = %lu\n",
		    n_violated_points, n_valid_points,
		    n_violated_points + n_valid_points);
	    }
	}

	// Normalize weights.
	printf("5\n");
	normalize_points(w, n_points);

	printf("6\n");
	// Resample points according to weights to construct policy.
	size_t *policy_idx = malloc(sizeof(size_t) * policy_size);
	kitagawa_resampling(w, n_points, policy_idx, policy_size);
	for (size_t i = 0; i < policy_size; ++i)
	{
	    size_t curr_point_idx = ndims * i;
	    size_t curr_policy_point_idx = ndims * policy_idx[i];

	    for (size_t d = 0; d < ndims; ++d)
	    {
		policy[curr_point_idx+d] = grid_points[curr_policy_point_idx+d];
	    }

	    // Normalize mixture weights P_LL, P_LH, and P_HL.
	    normalize_points(policy+curr_point_idx+3, 3);
	}

	printf("7\n");
	// Free memory that needs to be freed.
	free(policy_idx);
	free(w);
	free(grid_points);
	free(npoints_per_dim);
    }
    printf("build_sir_policy(): END\n");
}

