#ifndef __UTILS_H_
#define __UTILS_H_

#include <stddef.h>

double compute_E_in
    (const int* S_in, const double* h_effective_in, const double* h,
     int num_nodes);

size_t compute_npoints_in_grid(size_t *npoints_per_dim, size_t ndims);

double* make_deterministic_grid_points_nd
    (double *bounds, size_t *npoints_per_dim, size_t ndims);

double* make_uniform_grid_points_nd
    (double *bounds, size_t *npoints_per_dim, size_t ndims, size_t npoints);

// All points must be non-negative.
void normalize_points(double *points, const size_t n_points);

#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#endif // __UTILS_H_

