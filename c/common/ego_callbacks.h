#ifndef _EGO_CALLBACKS_H__
#define _EGO_CALLBACKS_H__

#include <stdlib.h>

// Function pointer to the callback that maximizes the acquisition function.
// Parameters: bounds_ptr, max_x_ptr, x_size
typedef int (*acq_maximizer_fn_t)(double *, double *, size_t);

// Function pointer to the callback that adds a datapoint to the GP.
// Parameters: x_ptr, x_size, y 
typedef int (*gp_add_data_fn_t)(double *, size_t, double);

// Function pointer to the callback that finds the coordinates (x, y) that
// maximizes y in the GP.
// Parameters: x_ptr, x_size, y_ptr
typedef int (*gp_maximizer_fn_t)(double *, size_t, double *);

// Function pointer to the callback that returns the GP posterior mean for the
// coordinates (x, y).
// Parameters: coords_ptr, coords_size
typedef double (*gp_posterior_mean_fn_t)(double *, size_t);


// Function pointer to the callback that finds the coordinates (x, y) which
// is a local optima of the UCB rule.
// Parameters: x_ptr, x_size, y_ptr
typedef int (*gp_maximizer_local_fn_t)(double *, size_t, double *);


typedef int (*sample_gp_posterior_mean_fn_t)(double *, size_t, double *, double *, size_t);


typedef int (*return_sample_fn_t)(int *, size_t);
#endif // _EGO_CALLBACKS_H__

