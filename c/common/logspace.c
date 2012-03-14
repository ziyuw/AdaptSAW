#include "logspace.h"
#include "utils.h"

#include <math.h>
#include <float.h>

double log_sum_exp(double log_x, double log_y)
{
    double log_max = max(log_x, log_y);

    log_x -= log_max;
    log_y -= log_max;

    double result = log(exp(log_x) + exp(log_y)) + log_max;
    return result;
}

double max_n(double *nums, const size_t n_nums)
{
    double result = -DBL_MAX;

    for (size_t i = 0; i < n_nums; ++i)
    {
        if (nums[i] > result)
        {
            result = nums[i];
        }
    }

    return result;
}

double log_sum_exp_nitems(double *log_nums, const size_t n_log_nums)
{
    double log_max = max_n(log_nums, n_log_nums);

    for (size_t i = 0; i < n_log_nums; ++i)
    {
        log_nums[i] -= log_max;
    }

    double result = 0.0;
    for (size_t i = 0; i < n_log_nums; ++i)
    {
        result += exp(log_nums[i]);
    }

    result = log(result) + log_max;
    return result;
}

