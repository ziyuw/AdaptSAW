#include "utils.h"
#include <mt19937ar.h>

#include <stdlib.h>

#include <math.h>

double compute_E_in
    (const int* S_in, const double* h_effective_in, const double* h,
     int num_nodes)
{
    double E_in = 0.0;

    for (int i = 0; i < num_nodes; ++i)
    {
        E_in -= S_in[i]*h[i];
    }
    for (int i = 0; i < num_nodes; ++i)
    {
        E_in -= S_in[i]*(h_effective_in[i] - h[i])/2.0;
    }

    return E_in;
}

size_t compute_npoints_in_grid(size_t *npoints_per_dim, size_t ndims)
{
    size_t npoints = 1;
    for (size_t i = 0; i < ndims; ++i) { npoints *= npoints_per_dim[i]; }

    return npoints;
}

double* make_deterministic_grid_points_nd
    (double *bounds, size_t *npoints_per_dim, size_t ndims)
{
    size_t npoints = compute_npoints_in_grid(npoints_per_dim, ndims);

    // Allocate an array to store the points.
    double *points = malloc(sizeof(double) * npoints * ndims);

    // Compute increments of each dimension
    double *increments = malloc(sizeof(double) * ndims);
    for (size_t i = 0; i < ndims; ++i)
    {
        increments[i] =
            (bounds[2*i + 1] - bounds[2*i]) / (double)(npoints_per_dim[i] - 1);
    }

    // Compute increment intervals of each dimension
    size_t *incr_intervals = malloc(sizeof(size_t) * ndims);
    for (size_t i = 0; i < ndims; ++i) { incr_intervals[i] = 1; }

    if (ndims > 1)
    {
        for (size_t i = 1; i < ndims; ++i)
        {
            incr_intervals[i] = incr_intervals[i-1] * npoints_per_dim[i-1];
        }
    }

    // Set up current position per dimension 
    size_t *pos_per_dim = malloc(sizeof(size_t) * ndims);
    for (size_t i = 0; i < ndims; ++i) { pos_per_dim[i] = 0; }

    // Construct grid points by looping over indices in points array, applying
    // increments only to dimensions that match the current index. The current
    // index modded by the number of unique points on each dimension multiplied
    // with other dimensions determines whether the dimension matches the
    // current index, to have the increment applied.
    for (size_t p = 0; p < npoints; ++p)
    {
        size_t param_idx = p*ndims;
        for (size_t i = 0; i < ndims; ++i)
        {
            points[param_idx + i] = bounds[2*i] + increments[i]*pos_per_dim[i];
        }

        ++pos_per_dim[0];
        pos_per_dim[0] = pos_per_dim[0] % npoints_per_dim[0];

        if (p == 0) { continue; }

        // Increment the positions per dimension
        size_t incr_idx = ndims-1;
        for (; incr_idx >= 1; --incr_idx)
        {
            if (p % incr_intervals[incr_idx] == 0) { break; }
        }

        for (size_t i = 1; i <= incr_idx; ++i)
        {
            ++pos_per_dim[i];
            pos_per_dim[i] = pos_per_dim[i] % npoints_per_dim[i];
        }
    }

    // Free memory that should be freed.
    free(pos_per_dim);
    free(incr_intervals);
    free(increments);

    return points;
}

double* make_uniform_grid_points_nd
    (double *bounds, size_t *npoints_per_dim, size_t ndims, size_t npoints)
{
    // Allocate an array to store the points.
    double *points = malloc(sizeof(double) * npoints * ndims);

    // Generate points.
    for (size_t p = 0; p < npoints; ++p)
    {
        double *curr_point = points + p*ndims;

        // Generate SAW_length_[min, max] that satisfy constraints.
        {
            double SAW_length_min =
                bounds[2] + genrand_real1() * (bounds[3] - bounds[2]);

            double SAW_length_max_lb = max(SAW_length_min, bounds[4]);
            double SAW_length_max =
                SAW_length_max_lb +
                genrand_real1() * (bounds[5] - SAW_length_max_lb);

            curr_point[1] = SAW_length_min;
            curr_point[2] = SAW_length_max;
        }

        // Generate all other values from the uniform distribution.
        curr_point[0] = bounds[0] + genrand_real1() * (bounds[1] - bounds[0]);

        for (size_t i = 3; i < ndims; ++i)
        {
            curr_point[i] =
                bounds[2*i] + genrand_real1() * (bounds[2*i+1] - bounds[2*i]);
        }
    }

    return points;
}

#include "logspace.h"
#include <math.h>
void normalize_points(double *points, const size_t n_points)
{
    for (size_t i = 0; i < n_points; ++i)
    {
        points[i] = log(points[i]);
    }

    double log_Z = log_sum_exp_nitems(points, n_points);

    for (size_t i = 0; i < n_points; ++i)
    {
        points[i] -= log_Z;
        points[i] = exp(points[i]);
    }
}

