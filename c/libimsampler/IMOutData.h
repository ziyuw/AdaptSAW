#ifndef __IMOUTDATA_H_
#define __IMOUTDATA_H_

struct _IMMCMCDiagnostics {

  int n_moves;
  double* E_samples;   // n_moves
  double* E_proposed;   // n_moves
  double* log_f_fwd; // n_moves
  double* log_f_rev; // n_moves
  double* MH_ratio;   // n_moves

  int *SAW_length; // n_moves
  double *gammas; // n_moves

  double *rewards; // n_moves

  int n_actions;

  int *S_end;
};
typedef struct _IMMCMCDiagnostics IMMCMCDiagnostics;

void allocMCMCDiagnostics(IMMCMCDiagnostics* p_data, int n_moves);
void trashMCMCDiagnostics(IMMCMCDiagnostics* p_data);
void storeDiagnosticsData
    (IMMCMCDiagnostics* p_data, int curr_move, double E_0, double E_prime,
     double log_f_fwd, double log_f_rev, double MH_ratio, int SAW_length,
     double gamma);

void storeEndState(IMMCMCDiagnostics* p_data, int *S_end, int num_nodes);

#endif // __IMOUTDATA_H_

