#ifndef __SARDOUTDATA_H_
#define __SARDOUTDATA_H_

#include "SARD.h"

struct _SARDMCMCDiagnostics {
  int n_moves;

  // 2D arrays all manually implemented as COLUMN-MAJOR 
  double* E_samples;   // n_moves
  double* E_proposed;   // n_moves
  double* log_f_fwd; // n_moves
  double* log_f_rev; // n_moves
  double* MH_ratio;   // n_moves
  int* move_type; // n_moves

  int *used_sigma_lengths; // n_moves x iters_per_move
  int *used_rho_lengths; // n_moves x iters_per_move
  int *iters_per_move; // n_moves

  double *used_gammas_his; // n_moves
  double *used_gammas_los; // n_moves

  double *used_P_LL; // n_moves
  double *used_P_LH; // n_moves
  double *used_P_HL; // n_moves

  double *rewards; // n_moves

  int *S_end; // numnodes

  int nactions;
};
typedef struct _SARDMCMCDiagnostics SARDMCMCDiagnostics;

void allocMCMCDiagnostics(SARDMCMCDiagnostics* p_data, int n_moves);
void trashMCMCDiagnostics(SARDMCMCDiagnostics* p_data);
void storeDiagnosticsData
    (SARDMCMCDiagnostics* p_data, const SARDSimParams *p_sim_params,
     int curr_move, double E_0, double E_prime,
     double log_f_fwd, double log_f_rev,
     double MH_ratio, int move_type,
     int *sigma_lengths, int *rho_lengths, int iters_per_move,
     double used_gamma_hi, double used_gamma_lo,
     double P_LL, double P_LH, double P_HL, double reward);

void storeEndState(SARDMCMCDiagnostics* p_data, int *S_end, int num_nodes);

#endif // __SARDOUTDATA_H_

