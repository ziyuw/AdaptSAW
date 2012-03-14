from ego.gaussianprocess import GaussianProcess
from ego.acquisition import maximizeEI
from ego.acquisition import maximizeUCB_local
from ego.acquisition import maximizeMean

from ctypes import *
from numpy import *

# Definitions of c function pointer types to the three functions
# needed to implement EGO:
#   - a function to maximize the acquisition function (maximizeEI).
#   - a function to add a new data point to the GP.
#   - a function to get the (x, y) with maximum y in the GP.

# ctypes function pointer for a function that maximizes an acquisition
# function.
ACQ_MAXIMIZER_FUNC = \
     CFUNCTYPE(# Nothing is returned, need to change this to void
               c_int,
               # Bounds
               POINTER(c_double),
               # Maximizing x-coordinates and coordinates size
               POINTER(c_double), c_size_t)

# I am not sure how to pass the GP into C through ctypes, so the most
# obvious thing to do is probably to create a closure over the GP,
# hold a reference to the closure and create a function pointer to the
# closure.
def maximizeEI_wrap(gp, bounds_ptr, max_x_ptr, x_size):
    bounds = [[bounds_ptr[2*i], bounds_ptr[2*i+1]] for i in xrange(x_size)]
    #print bounds
    _, max_x = maximizeEI(gp, bounds)

    assert(len(max_x) == x_size)

    #print max_x
    for i in xrange(x_size):
        max_x_ptr[i] = max_x[i]
    x_size = c_size_t(len(max_x))

    return 0

# ctypes function pointer for a function that adds a new data point to
# the GP.
GP_ADD_DATA_FUNC = \
     CFUNCTYPE(# Nothing is returned, need to change this to void
               c_int,
               # x-coordinates, coordinates size and y
               POINTER(c_double), c_size_t, c_double)

# This one has to be wrapped in a closure with the GP as well.
def add_data_to_gp(gp, x_ptr, x_size, y):
    #print x_ptr[:x_size], y
    gp.addData(x_ptr[:x_size], y)
    return 0

# ctypes function pointer for a function that gets the (x, y) with
# maximum y in the GP
GP_MAXIMIZER_FUNC = \
    CFUNCTYPE(# Nothing is returned, need to change this to void
              c_int,
              # Output parameters x and y
              POINTER(c_double), c_size_t, POINTER(c_double))

def maximize_gp(gp, x_ptr, x_size, y_ptr):
    max_idx = argmax(gp.Y)
    x, y = gp.X[max_idx], gp.Y[max_idx]

    assert(len(x) == x_size)

    for i in xrange(x_size):
        x_ptr[i] = x[i]
    y_ptr[0] = y
    
    return 0

GP_POSTERIOR_MEAN_FUNC = \
    CFUNCTYPE(c_double, # value of gp posterior mean
              POINTER(c_double), # coordinates
              c_size_t) # size of coordinates

def gp_posterior_mean(gp, coords_ptr, coords_size):
    coords = zeros(coords_size)
    for i in xrange(coords_size):
        coords[i] = coords_ptr[i]
        
    return gp.posterior(coords)[0]

LOCAL_ACQ_MAXIMIZER_FUNC = \
     CFUNCTYPE(# Nothing is returned, need to change this to void
               c_int,
               # Bounds
               POINTER(c_double),
               # Maximizing x-coordinates and coordinates size
               POINTER(c_double), c_size_t)


def maximizeEI_local_wrap(gp, bounds_ptr, max_x_ptr, x_size, strategy):  
    """
	method can either be set to 'direct' or 'l-bfgs-b'
    """    
    if strategy == 1 or strategy == 4:
	method = 'l-bfgs-b'
    elif strategy == 3 or strategy == 2:
	method = 'direct'

    # Put bounds in the format required by l-bfgs-b
    bounds = [(bounds_ptr[2*i], bounds_ptr[2*i+1]) for i in xrange(x_size)]
    
    # Initialize x0
    import random as ran
    
    x0 = zeros(x_size)    
    
    for i in xrange(x_size):
	x0[i] = ran.uniform(bounds[i][0], bounds[i][1])

    # Maximize the function
    if method == 'direct':
	max_x, max_val = maximizeMean(gp, bounds)
    elif method == 'l-bfgs-b':
	max_x, max_val = maximizeUCB_local(gp, bounds, x0)

    assert(len(max_x) == x_size)

    for i in xrange(x_size):
        max_x_ptr[i] = max_x[i]
    
    # Store the objective value at the end
    max_x_ptr[x_size] = -max_val
    
    x_size = c_size_t(len(max_x)+1)

    return 0

SAMPLE_GP_POSTERIOR_MEAN_FUNC = \
    CFUNCTYPE(c_double, # value of gp posterior mean
              POINTER(c_double), # coordinates
              c_size_t, # size of the coordinate
              POINTER(c_double), # bounds
              POINTER(c_double), # pointer to the samples acquired
              c_size_t # Number of samples needed
              ) 
              
def sample_gp_posterior_mean(gp, coords_ptr, coords_size, bounds_ptr, \
			     sample_x_ptr, num_samples, check_constraint = None, exp_temp = 10, MH_rate = 0.09):
				 
    """
	Sample by using the Metropolis method to get policies in the re-sampling phase.
	We treat the exponentiated (GP-mean*exp_temp) as unnormalized density.
	The proposal distribution is normal distribution: N(Mu, MH_rate*I).
	
    """
    bounds = [(bounds_ptr[2*i], bounds_ptr[2*i+1]) for i in xrange(coords_size)]
    
    def check_violate_constraint(nx_pt, bounds, coords_size, check_constraint = None):
	for i in xrange(coords_size):
	    if not (nx_pt[i] <= bounds[i][1] and nx_pt[i] >= bounds[i][0]):
		return True
	if not check_constraint is None:
	    return check_constraint(nx_pt, bounds, coords_size)
	return False


    coords = zeros(coords_size)
    for i in xrange(coords_size):
        coords[i] = coords_ptr[i]
    
    import numpy.random as nr
    
    cur_pt = coords
    
    factor = MH_rate
    cov = eye(coords_size)
    for i in xrange(coords_size):
	cov[i,i] = (bounds[i][1] - bounds[i][0])**2*factor;
    
    for i in xrange(num_samples):
	for j in xrange(coords_size):
	    sample_x_ptr[i*coords_size+j] = cur_pt[j]
	
	nx_pt = nr.multivariate_normal(cur_pt, cov)
	while check_violate_constraint(nx_pt, bounds, coords_size, check_constraint):
	    nx_pt = nr.multivariate_normal(cur_pt, cov)

	mh_ratio = exp((- gp.posterior(cur_pt)[0] + gp.posterior(nx_pt)[0])*exp_temp)
	r = nr.uniform()
	
	if r < mh_ratio:
	    cur_pt = nx_pt
	        
    return 0
    
RETURN_SAMPLE = \
    CFUNCTYPE(c_double, # value of return value
              POINTER(c_int), # coordinates
              c_size_t # Dimenion of sample
              )
              
def return_sample(samples, sample, dim):
    py_sample = zeros(dim)
    
    
    for i in xrange(dim):
	py_sample[i] = sample[i]
    
    samples.append(py_sample)
    
    return 0
    