from __future__ import division
import tests
import ising as isg

import sys 
from numpy import *
import matplotlib.pyplot as plt
import utils

"""
Example Usage:

Runs SARDONICS on the ferro meganetic model with adaptation
    python run_sard_tests.py ferro 

Runs SARDONICS on the ferro meganetic model without adaptation
    python run_sard_tests.py ferro noadapt
"""

noadapt = False
if len(sys.argv) > 2:
    test_name = sys.argv[1]
    noadapt = bool(sys.argv[2])
elif len(sys.argv) > 1:
    test_name = sys.argv[1]
else:
    test_name = 'ferro'

# Run tests based on the tests in the UAI 2010 paper.
if test_name == 'ferro':
    print "Running Ferro2D_UAI."
    test = isg.Ferro2D_UAI(60, 60)
elif test_name == 'isg3d':
    print "Running ISG3D_UAI."
    test = isg.ISG3D_UAI(9, 9, 9)
elif test_name == 'rbm':
    print "Running RBM_UAI."
    test = isg.RBM_UAI()
elif test_name == 'ferrof':
    print "Running Ferro2D_Frustrated."
    test = isg.Ferro2D_Frustrated(60, 60)
elif test_name == 'chi':
    print "Running Chimera."
    test = isg.test = isg.Chimera(4,4,4)

if noadapt:
    #=====================================================
    # Run SARDONICS without parameter adaptation
    #=====================================================
    
    # strategy: 0 choose parameters uniformly within bounds
    #	        1 choose parameters deterministicly 
    
    testspec_sard_unif = \
	{"graph" : test,
	"nmoves" : 8*10**4,
	"gamma_hi_bounds" : array([1.0, 2.5]) * test.beta_true,
	"gamma_lo_bounds" : array([1.0, 1.5]) * test.beta_true,
	"max_SAW_length_bounds" : [60, 61], "min_SAW_length_bounds" : [50, 51],
	"iters_per_move_bounds" : [1, 2], "strategy" : 1, 'P_LL' : 1.0, 'P_HL' : 0.0,}

    sard_unif_test = tests.SARD_Unif_Test(**testspec_sard_unif)
    sard_unif_test.run()

    plt.figure()
    plt.subplot(311)
    line1 = plt.plot(sard_unif_test.MH_ratios, 'b*')

    plt.subplot(312)
    line2 = plt.plot(sard_unif_test.E_samples, '-b')
    plt.grid(True)
    
    
    max_lag = 500
    plt.subplot(313)
    plt.plot(map(lambda i: utils.compute_acf(sard_unif_test.E_samples[1:], i), xrange(max_lag+1)))
    plt.grid(True)
    plt.ylim(ymax=1.0)
    plt.ylim(ymin=0.0)
    
    plt.show()
    
else:
    #=====================================================
    # Run SARDONICS with parameter adaptation
    #=====================================================
    
    # Total number of moves (including the moves during adaptation)
    nmoves = 10**5

    # Set total number of adaptation
    nadaptations = 200
    
    # Set the number of moves used to calculate the reward for each adaptation step
    nruns_per_param_update = 1*10**2
    
    # Set the policy size
    policy_size = 1000
    
    """
    'strategy' specifying what to do for the "resampling phase"
		0 Resampling through deterministic grid points (Only works with low dimensions and not recommended)
    	1 multiple start l-bfgs-b + Metropolis sampling of the GP posetior
    	2 direct + Metropolis sampling of the GP posetior
    """
    strategy = 2
    
    """
    'return_samples' specifies whether samples are to be returned
    	0: do not return
    	1: return samples after adaptation is finished
    	2: return all samples
    """
    return_samples = 0
    
    # Set the bounds for Self-Avoiding-Walk (SAW) length
    min_SAW_length_bounds = [1, 70]
    max_SAW_length_bounds = [2, 120]
    
    # Set the bounds for number of iterated SAW
    iters_per_move_bounds = [1, 5]
    
    # Set the bounds for gammas 
    gamma_lo_bounds  = array([0.89, 1.05]) * test.beta_true
    gamma_hi_bounds = array([0.9, 1.15]) * test.beta_true
    
    
    testspec_sard_ego = \
	{"graph" : test,
	"nmoves" : nmoves, "nadaptations" : nadaptations,
	"nruns_per_param_update" : nruns_per_param_update, "policy_size" : policy_size,
	"gamma_hi_bounds" : gamma_hi_bounds, "gamma_lo_bounds" : gamma_lo_bounds,
	"max_SAW_length_bounds" : max_SAW_length_bounds, "min_SAW_length_bounds" : min_SAW_length_bounds,
	"iters_per_move_bounds" : iters_per_move_bounds, "strategy" : strategy, "exp_temp" : 10, 
	"MH_rate" : 0.005**2, "return_samples" : return_samples,}

    testspec_sard_ego_save = {}
    for key in testspec_sard_ego:
	if not key == 'graph':
	    testspec_sard_ego_save['run_'+key] = testspec_sard_ego[key]

    # Set SARD parameters
    sard_ego_test = tests.SARD_EGO_Test(**testspec_sard_ego)
    
    import time
    time_1 = time.time()
    
    # Run SARDONICS with adaptation
    sard_ego_test.run()
    time_2 = time.time()
    print "Sampling lasted: ", time_2 - time_1
    
    # To get the samples call 'sard_ego_test.samples'
    
    # Save results and save plots
    sard_ego_test.save_results(other_specs = testspec_sard_ego_save)
    sard_ego_test.plot_acf(start = nruns_per_param_update*nadaptations)
    sard_ego_test.plot_E_samples()

    # Plot the results if show_plot is set to be True
    show_plot = False
    if show_plot:
	plt.figure(1)
	plt.subplot(221)
	line1 = plt.plot(sard_ego_test.MH_ratios, 'b*')

	plt.subplot(223)
	line2 = plt.plot(sard_ego_test.E_samples, '-b*')
	plt.grid(True)

	plt.subplot(222)
	line3 = plt.plot(sard_ego_test.rewards, '-r.')

	plt.subplot(224)
	max_lag = 500
	plt.plot(map(lambda i: utils.compute_acf(sard_ego_test.E_samples[nruns_per_param_update*nadaptations:], i), xrange(max_lag+1)))
	plt.ylim(ymax=1.0)
	plt.ylim(ymin=0.0)

	plt.grid(True)
	plt.show()
