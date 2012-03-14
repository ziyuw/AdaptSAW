import tests

nmoves = 20000
nadaptations = 100
nruns_per_param_update = 100

nadapted_moves = tests.get_num_adapted_moves(nruns_per_param_update, nadaptations)

def plot_EGO_based_test(testspec, testnum=None):
    testtype, args = testspec
    test = testtype(*args)
    nresults = test.get_num_results()

    if testnum is None:
        testnum = 0

    test.load_results(testnum)

    test.plot_gp_mean()
    test.plot_gp_query_points()
    test.plot_mh_ratio_moving_avg()

import matplotlib.pyplot as plt
import re
def plot_energies(testspecs, plot_fn, testnum=None):
    if testnum is None:
        testnum = 0
        
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for testspec in testspecs:
        testtype, args = testspec
        test = testtype(*args)
        test.load_results(testnum)

        E_samples = test.E_samples
        if re.match('IM_EGO', test.name):
            E_samples = E_samples[10000:14000]
        else:
            E_samples = E_samples[:4000]

        test_name = re.sub('_', '', test.name)
        test_name = re.sub('EGO', 'BayesOpt', test_name)
        ax.plot(E_samples, label=test_name)

    ax.legend()
    ax.set_xlabel('Steps')
    ax.set_ylabel('Energy')
    
    fig.savefig(plot_fn, format='pdf')

from numpy import *

def get_energy(testtype, args, count):
    test = testtype(*args)
    res_num = test.get_last_res_num() - count
    test.load_results(res_num)

    E_samples = test.E_samples
    if re.match('IM_EGO', test.name):
	E_samples = E_samples[10000:14000]
    else:
	E_samples = E_samples[:4000]
	
    return E_samples
    

def plot_avg_energy(test):
    
    num_samplers = 5
    
    total_energy = zeros(4000)
    
    for count in xrange(1, num_samplers+1):
	energy = get_energy(test[0], test[1], count).T
	total_energy = total_energy + energy
	
    average_energy = total_energy/num_samplers
    
    testtype = test[0]
    args = test[1]
    test = testtype(*args)
    
    print average_energy, average_energy.shape
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    test_name = re.sub('_', '', test.name)
    test_name = re.sub('EGO', 'BayesOpt', test_name)
    ax.plot(average_energy[0], label=test_name)

    print 'get here?'

    ax.legend()
    ax.set_xlabel('Steps')
    ax.set_ylabel('Energy')
    
    fig.savefig('something', format='pdf')




















