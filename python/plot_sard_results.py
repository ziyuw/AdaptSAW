from __future__ import division
import tests
import ising as isg

import sys 
from numpy import *
import matplotlib.pyplot as plt
import utils


unif = False
if len(sys.argv) > 2:
    test_name = sys.argv[1]
    unif = bool(sys.argv[2])
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


if unif:
    testspec_sard_unif = \
	{"graph" : ferro2d,
	"nmoves" : 10**3,
	"gamma_hi_bounds" : array([1.5, 2.5]) * ferro2d.beta_true,
	"gamma_lo_bounds" : array([0.5, 1.5]) * ferro2d.beta_true,
	"max_SAW_length_bounds" : [5, 15], "min_SAW_length_bounds" : [1, 10],
	"iters_per_move_bounds" : [1, 10]}

    sard_unif_test = tests.SARD_Unif_Test(**testspec_sard_unif)
    sard_unif_test.run()

    plt.figure(2)
    plt.subplot(211)
    line1 = plt.plot(sard_unif_test.MH_ratios, 'b*')

    plt.subplot(212)
    line2 = plt.plot(sard_unif_test.E_samples, '-b*')
    plt.grid(True)

    plt.grid(True)
    plt.show()
    
else:
    max_SAW_length_bounds_orig = [6, 15]
    min_SAW_length_bounds_orig = [1, 10]
    iters_per_move_bounds_orig = [1, 10]
    gamma_hi_bounds_orig = array([1.5, 2.5]) * test.beta_true
    gamma_lo_bounds_orig  = array([0.5, 1.5]) * test.beta_true
    
    nmoves = 10**5
    nadaptations = 200
    nruns_per_param_update = 1*10**2
    policy_size = 1000
    max_SAW_length_bounds = [20, 50]
    min_SAW_length_bounds = [1, 15]
    iters_per_move_bounds = [1, 10]
    gamma_hi_bounds = array([1.5, 2.5]) * test.beta_true
    gamma_lo_bounds  = array([0.5, 1.5]) * test.beta_true
    
    testspec_sard_ego = \
	{"graph" : test,
	"nmoves" : nmoves, "nadaptations" : nadaptations,
	"nruns_per_param_update" : nruns_per_param_update, "policy_size" : policy_size,
	"gamma_hi_bounds" : gamma_hi_bounds, "gamma_lo_bounds" : gamma_lo_bounds,
	"max_SAW_length_bounds" : max_SAW_length_bounds, "min_SAW_length_bounds" : min_SAW_length_bounds,
	"iters_per_move_bounds" : iters_per_move_bounds, "strategy" : 2, "exp_temp" : 10, "MH_rate" : 0.005**2,}

    testspec_sard_ego_save = {}
    for key in testspec_sard_ego:
	if not key == 'graph':
	    testspec_sard_ego_save['run_'+key] = testspec_sard_ego[key]

    sard_ego_test = tests.SARD_EGO_Test(**testspec_sard_ego)
    
    
    sard_ego_test.load_results()
    #sard_ego_test.plot_acf(start = nruns_per_param_update*nadaptations)
    #sard_ego_test.plot_E_samples()
    sard_ego_test.plot_policy_as_histogram()

    #plt.figure(1)
    #plt.subplot(221)
    #line1 = plt.plot(sard_ego_test.MH_ratios, 'b*')

    #plt.subplot(223)
    #line2 = plt.plot(sard_ego_test.E_samples, '-b*')
    #plt.grid(True)

    #plt.subplot(222)
    #line3 = plt.plot(sard_ego_test.rewards, '-r.')

    #plt.subplot(224)
    #max_lag = 500
    #plt.plot(map(lambda i: utils.compute_acf(sard_ego_test.E_samples[nruns_per_param_update*nadaptations:], i), xrange(max_lag+1)))
    #plt.ylim(ymax=1.0)
    #plt.ylim(ymin=0.0)

    #plt.grid(True)
    #plt.show()
