#ifndef __SARDUTILS_H_
#define __SARDUTILS_H_

void initSARDState(SARDState* p_sard_state,
                   const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params );
void setSARDState(SARDState* p_sard_state,
                const SARDInputModel* p_in_model,
                const int *S_new, const double* h_effective_new);
void freeSARDState(SARDState* p_sard_state,
                   const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params );
void initSAWData(SAWData* p_saw_data,
                 const SARDInputModel* p_in_model,
                 const SARDSimParams* p_sim_params);
void freeSAWData(SAWData* p_saw_data,
                 const SARDInputModel* p_in_model,
                 const SARDSimParams* p_sim_params );

void resetSAWHeaps(SAWData* p_saw_data, const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params);

void resetSAWHeapsCurrent(SARDState* p_sard_state, SAWData* p_saw_data, const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params);

void proposeSAWMove(SAWData* p_saw_data,
                    SARDState* p_sard_state,
                    const SARDInputModel* p_in_model,
                    const SARDSimParams* p_sim_params);
void undoSAWMove(SAWData* p_saw_data,
                 SARDState* p_sard_state,
                 const SARDInputModel* p_in_model,
                 const SARDSimParams* p_sim_params);

void evalLogFs(double* p_log_f_fwd, double* p_log_f_rev,
               const SAWData* p_saw_data, const SARDSimParams* p_sim_params);

#endif // __SARDUTILS_H_

