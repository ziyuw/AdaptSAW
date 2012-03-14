#include "objectives.h"
#include "utils.h"

#include <math.h>
#include <float.h>

double compute_mean(double *numbers, int n)
{
    double mean = 0.0;

    for (int i = 0; i < n; ++i)
    {
        mean += numbers[i];
    }

    mean /= (double)n;
    return mean;
}

double compute_variance(double *numbers, int n, double mean)
{
    double var = 0;

    for (int i = 0; i < n; ++i)
    {
        double x = numbers[i] - mean;
        x *= x;
        var += x;
    }

    var /= (double)n;
    return var;
}

// Computes the auto-correlation over n samples with a k-step lag, starting
// from E_samples[0] and going up to E_samples[n - k - 1].
double compute_acf(double *E_samples, int n, int k,
                   double mean, double variance)
{
    if (variance == 0.0)
    {
        variance = DBL_MIN;
    }

    double acf = 0.0;

    for (int t = 0; t < n - k; ++t)
    {
        double x = (E_samples[t] - mean) * (E_samples[t+k] - mean);

        if (x == 0.0)
        {
            x = DBL_MIN;
        }

        acf += x;
    }

    acf /= variance;
    acf /= (double)(n - k);

    // This can happen due to rounding errors.
    if (acf > 1)
    {
        acf = 1;
    }
    else if (acf < -1)
    {
        acf = -1;
    }

    return acf;
}

// Computes the average auto-correlation over all possible step-lengths, up to
// a maximum step-length of 500.
double compute_avg_abs_acf(double *E_samples, int nsamples, int max_lag_size)
{
    double mean = compute_mean(E_samples, nsamples);
    double var = compute_variance(E_samples, nsamples, mean);

    double avg_abs_acf = 0.0;

    int min_k = 1;
    int max_k = min(max_lag_size, min(500, nsamples - 1));

    for (int k = min_k; k <= max_k; ++k)
    {
        avg_abs_acf += fabs(compute_acf(E_samples, nsamples, k, mean, var));
    }

    avg_abs_acf /= (double)(max_k - min_k + 1);
    return avg_abs_acf;
}

double compute_robust_avg_abs_acf_est
    (double *E_samples,
     int nruns_per_param_update, int first_run_with_curr_param_idx)
{
    double reward = 0.0;

    for (int i = 0; i < nruns_per_param_update-25; ++i)
    {
        double *acf_window_start =
            E_samples+first_run_with_curr_param_idx+i;
        int acf_window_size = nruns_per_param_update-i;

        double avg_abs_acf =
            compute_avg_abs_acf(acf_window_start, acf_window_size, 100);
        reward += 1.0 - avg_abs_acf;
    }

    reward /= nruns_per_param_update-25;

    return reward;
}

