#ifndef __OBJECTIVES_H_
#define __OBJECTIVES_H_

#include <stddef.h>

typedef enum objectives
{
    MH_RATIO_AVG_ABS_DIFF,
    MH_RATIO,
    ACF
}
objective;

#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
 
double compute_mean(double *numbers, int n);
double compute_variance(double *numbers, int n, double mean);
double compute_acf(double *E_samples, int nsamples,
                   int k, double mean, double variance);

double compute_avg_abs_acf(double *E_samples, int nsamples, int max_lag_size);

double compute_robust_avg_abs_acf_est
    (double *E_samples,
     int nruns_per_param_update, int first_run_with_curr_param_idx);

#endif // __OBJECTIVES_H_

