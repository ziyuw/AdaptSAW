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
#include "ego_callbacks.h"

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

void sampleSARD_unif(const SARDInputModel* p_in_model,
                     SARDSimParams* p_sim_params,
                     SARDMCMCDiagnostics* p_out_data,
		     return_sample_fn_t return_sample_fn,
		     double P_LL, double P_HL, size_t strategy)
{
    SAWData* p_saw_data = (SAWData*)malloc(sizeof(SAWData));
    SARDState* p_sard_state = (SARDState*)malloc(sizeof(SARDState));

    initSARDState(p_sard_state, p_in_model, p_sim_params);
    // Set up the heaps, etc. 
    initSAWData(p_saw_data, p_in_model, p_sim_params);

    int nruns_until_next_param_update = 0;

    if (strategy == 1) {
	p_sim_params->SAW_length_min = \
	    p_sim_params->min_SAW_length_bound_min ;
	
	p_sim_params->SAW_length_max = p_sim_params->max_SAW_length_bound_min;

	p_sim_params->P_LL = P_LL;
	p_sim_params->P_LH = 1.0 - P_LL- P_HL;
	p_sim_params->P_HL = P_HL;

	// Pick iters_per_move, gamma_hi and gamma_lo.
	p_sim_params->iters_per_move =
	    p_sim_params->iters_per_move_min;
	p_sim_params->gamma_hi =
	    p_sim_params->gamma_hi_min;
	p_sim_params->gamma_lo =
	    p_sim_params->gamma_lo_min;
	
	    
	int debug = 1;
	if (debug) {
	    printf("iters_per_move: %d\n", (int)round(p_sim_params->iters_per_move));
	    printf("length_min: %d\n", p_sim_params->SAW_length_min);
	    printf("length_max: %d\n", p_sim_params->SAW_length_max);
	    
	    printf("p_sim_params->P_LL: %f\n", p_sim_params->P_LL);
	    printf("p_sim_params->P_LH: %f\n", p_sim_params->P_LH);
	    printf("p_sim_params->P_HL: %f\n", p_sim_params->P_HL);
	    
	    printf("p_sim_params->gamma_hi: %f\n", p_sim_params->gamma_hi);
	    printf("p_sim_params->gamma_lo: %f\n", p_sim_params->gamma_lo);
	}
	resetSAWHeaps(p_saw_data, p_in_model, p_sim_params);
    }
    
    for (int curr_move=0; curr_move<p_sim_params->n_moves; curr_move++)
    {
        if (nruns_until_next_param_update <= 0 && strategy == 0)
        {
            nruns_until_next_param_update = p_sim_params->nruns_per_param_update;
	    	    
            // Pick SAW_length_min and SAW_length_max, satisfying the
            // constraint that SAW_length_max >= SAW_length_min.
            p_sim_params->SAW_length_min = \
                p_sim_params->min_SAW_length_bound_min +
                (int) (genrand_real2() *
                       (p_sim_params->min_SAW_length_bound_max -
                        p_sim_params->min_SAW_length_bound_min + 1));
            int SAW_length_max_lb = \
                max(p_sim_params->SAW_length_min,
                    p_sim_params->max_SAW_length_bound_min);
            p_sim_params->SAW_length_max = \
                SAW_length_max_lb +
                (int) (genrand_real2() *
                       (p_sim_params->max_SAW_length_bound_max -
                        SAW_length_max_lb + 1));

            // Pick P_LL, P_LH and P_HL, satisfying the constraint that
            // P_LL + P_LH + P_HL = 1 by normalizing.
            p_sim_params->P_LL = genrand_real1();
            p_sim_params->P_LH = genrand_real1();
            p_sim_params->P_HL = genrand_real1();

            double mix_Z = \
                exp(log_sum_exp(log_sum_exp(log(p_sim_params->P_LL),
                                            log(p_sim_params->P_HL)),
                                log(p_sim_params->P_LH)));
            p_sim_params->P_LL /= mix_Z;
            p_sim_params->P_LH /= mix_Z;
            p_sim_params->P_HL /= mix_Z;

            // Pick iters_per_move, gamma_hi and gamma_lo.
            p_sim_params->iters_per_move =
                p_sim_params->iters_per_move_min +
                (int) (genrand_real2() *
                       (p_sim_params->iters_per_move_max -
                        p_sim_params->iters_per_move_min + 1));
            p_sim_params->gamma_hi =
                p_sim_params->gamma_hi_min +
                (genrand_real1() *
                 (p_sim_params->gamma_hi_max - p_sim_params->gamma_hi_min));
            p_sim_params->gamma_lo =
                p_sim_params->gamma_lo_min +
                (genrand_real1() *
                 (p_sim_params->gamma_lo_max - p_sim_params->gamma_lo_min));
	    resetSAWHeaps(p_saw_data, p_in_model, p_sim_params);
        }
	
        --nruns_until_next_param_update;

        double log_f_fwd, log_f_rev, MH_ratio;

        double E_initial= p_sard_state->E;

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
        
        
//      printf("Just before SAW_move\n");
        proposeSAWMove(p_saw_data, p_sard_state, p_in_model, p_sim_params);
// 	printf("Just after SAW_move\n");
        evalLogFs(&log_f_fwd, &log_f_rev, p_saw_data, p_sim_params);
	
        MH_ratio = p_sim_params->beta_true*(E_initial - p_sard_state->E);
        MH_ratio += (log_f_rev - log_f_fwd);
        MH_ratio = exp(MH_ratio);
        MH_ratio = MH_ratio >1.0 ? 1.0 : MH_ratio;

	
        double reward = 0.0;

	
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
        
	
        return_sample_fn(p_sard_state->S, p_in_model->num_nodes);
        
        if ((curr_move+1)%1000 == 0) {
	    printf("Iteration: %d has finished.\n", curr_move+1);
	}
        
    } // end moves loop

    storeEndState(p_out_data, p_sard_state->S, p_in_model->num_nodes);

    freeSARDState(p_sard_state, p_in_model, p_sim_params);
    freeSAWData(p_saw_data, p_in_model, p_sim_params);  

    free(p_sard_state);
    free(p_saw_data);
}

