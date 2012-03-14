#ifndef __SARD_UNIF_H_
#define __SARD_UNIF_H_

#include "ego_callbacks.h"
#include "SARD.h"

void sampleSARD_unif(const SARDInputModel* p_in_model,
                     SARDSimParams* p_sim_params,
                     SARDMCMCDiagnostics* p_out_data,
		     return_sample_fn_t return_sample_fn,
		     double P_LL, double P_HL, size_t strategy);

#endif // __SARD_UNIF_H_

