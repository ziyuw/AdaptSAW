#include "sampling.h"
#include <mt19937ar.h>
#include <math.h>
#include <float.h>
#include "logspace.h"

#include <stdio.h>

void build_multinomial_cdf(double *cdf, double *pdf, size_t pdf_size)
{
    double Z = 0;
    for (size_t i = 0; i < pdf_size; ++i)
    {
        Z += pdf[i];
    }

    cdf[0] = pdf[0] / Z;

    for (size_t i = 1; i < pdf_size; ++i)
    {
        cdf[i] = cdf[i-1] + pdf[i] / Z;
    }

    cdf[pdf_size - 1] = 1.0; // could be ~1.0 due to fp errors
}

size_t sample_multinomial_cdf(double *cdf, size_t cdf_size)
{
    double x = genrand_real2();

    // linear search to find the corresponding idx
    for (size_t i = 0; i < cdf_size; ++i)
    {
        if (x <= cdf[i])
        {
            return i;
        }
    }

    printf("ERR: sample_multinomial_cdf() failed\n");
    return -1;
}

size_t sample_uniform(size_t pdf_size)
{
    double x = genrand_real2();
    size_t result = (size_t)floor(x * (double)pdf_size);

    if (result == pdf_size)
    {
        result = pdf_size - 1;
    }

    return result;
}

// A stupid, basic and crappy linear search-based sampling algorithm.
int sample_from_discrete_log_pd(double *log_pd, int log_pd_size)
{
    double log_x = log(genrand_real2());
    double log_F = log(DBL_MIN);
    int i = 0;

    int sampled_idx = -1;

    for (i = 0; i < log_pd_size; ++i)
    {
        log_F = log_sum_exp(log_F, log_pd[i]);

        if (log_F >= log_x)
        {
            sampled_idx = i;
            break;
        }
    }

    if (sampled_idx < 0)
    {
        if (exp(log_F) < 1)
        {
            // This is due to rounding errors resulting in log_pd not
            // log-exp-summing up to exactly 1.0, so allocate that round-off
            // error to the last element in the multinomial pd.
            return log_pd_size - 1;
        }
        else
        {
            // OMG WTF this is badness.
            // Should really be throwing an exception, but whatevs.
            return -1;
        }
    }

    return sampled_idx;
}

void kitagawa_resampling
    (double *w, size_t n_weights, size_t *s_res, size_t n_res)
{
    // deterministic Kitagawa resampling
    double *u = malloc(sizeof(double) * n_res);
    double x = genrand_real2();
    for (size_t i = 0; i < n_res; ++i)
    {
        u[i] = ((double)i + x) / (double)(n_res);
    }

    double *cdf = malloc(sizeof(double) * n_weights);
    build_multinomial_cdf(cdf, w, n_weights);
    for (size_t i = 0; i < n_res; ++i)
    {
        size_t idx = n_res;
        for (size_t j = 0; j < n_weights; ++j)
        {
            if (u[i] <= cdf[j])
            {
                idx = j;
                break;
            }
        }

        s_res[i] = idx;
    }

    free(u);
    free(cdf);
}

