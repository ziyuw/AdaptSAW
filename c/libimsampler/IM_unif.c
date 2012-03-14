#include "IM.h"
#include "IMUtils.h"
#include <mt19937ar.h>
#include <math.h>
#include "IM_unif.h"

void 
sampleIM_unif(const IMInputModel * p_in_model,
              const IMSimParams * p_sim_params,
              IMMCMCDiagnostics * p_out_data)
{
    SAWData         saw_data;
    IMState         im_state;

    initIMState(&im_state, p_in_model, p_sim_params);

    double gamma_min = p_sim_params->gammas[0];
    double gamma_max = p_sim_params->gammas[p_sim_params->n_gammas - 1];

    // Set up the heaps, etc.
	initSAWData(&saw_data, p_in_model, p_sim_params);

    int nruns_until_next_param_update = 0;
    int n_accepted = 0;
    int n_proposed = 0;

    for (int curr_move = 0; curr_move < p_sim_params->n_moves; curr_move++)
    {
        if (nruns_until_next_param_update == 0)
        {
            nruns_until_next_param_update = p_sim_params->nruns_per_param_update;
            saw_data.SAW_length = p_sim_params->SAW_length_min;

            // Draw SAW length from NAIVE uniform distribution
            // U[SAW_length_min, SAW_length_max]
            saw_data.SAW_length =
                p_sim_params->SAW_length_min +
                (int) (genrand_real2() *
                       (p_sim_params->SAW_length_max -
                        p_sim_params->SAW_length_min + 1));

            // Draw gamma from NAIVE uniform distribution
            // U[gammas[0], gammas[n_gammas - 1]]
            saw_data.gamma_up =
                gamma_min +
                (genrand_real2() * (gamma_max - gamma_min));
            saw_data.gamma_down = saw_data.gamma_up;
        }

        --nruns_until_next_param_update;

        double E_initial = im_state.E;
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
}

