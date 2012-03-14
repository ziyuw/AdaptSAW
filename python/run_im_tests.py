from __future__ import division
import tests
import ising as isg

# TODO: (Code release for IM sampler)
#	1. Return data
#	2. Installation -- Later
#	3. Graph structure and parameter passing
#	4. Inerface -- First
#	5. Pass in values like MH_rate and exp_temp for sampling from GP posetior means -- DONE
#	6. Allow IM_unif and IM_Expert
#	7. Configuration file has to be explained
#	8. hdf5 file content has to be figured out!
#	9. Files for my use has to be removed!
#	10. Finish README.txt
#	11. What to do about Eric Brochu and Nim's copyright stuff

nadaptations = 100
nruns_per_param_update = 100

# Run tests based on the tests in the UAI 2010 paper.
#ferro2d = isg.Ferro2D_UAI(60, 60)
isg3d = isg.ISG3D_UAI(9, 9, 9)
#rbm = isg.RBM_UAI()

def create_tests_for_ising_model\
        (ising_model, expert_saw_lens, other_saw_lens,
         policy_size, nadaptations, nruns_per_param_update, 
         nsamples, strategy, exp_temp, MH_rate):
    kawasaki_params = (ising_model, nsamples)
    
    expert_params = list(kawasaki_params)
    expert_params.extend(expert_saw_lens)

    unif_params = list(kawasaki_params)
    unif_params.extend(other_saw_lens)

    ego_params = list(kawasaki_params)
    ego_params.extend([nadaptations, nruns_per_param_update])
    ego_params.extend([policy_size])
    ego_params.extend(other_saw_lens)
    ego_params.extend([strategy, exp_temp, MH_rate])
    
    tests_for_model = {
        'kawasaki' : [tests.Kawasaki_Sampler_Test, kawasaki_params],
        'im_unif' : [tests.IM_Unif_Test, unif_params],
        'im_expert' : [tests.IM_Expert_Test, expert_params],
        'im_ego' : [tests.IM_EGO_Test, ego_params],
        'im_ego_premixed' : [tests.IM_EGO_premixed_Test, ego_params]
        }
    
    return tests_for_model
"""
    Strategy specifying what to do for the "resampling phase"
    strategy:	0 Nim's Resampling
		1 l-bfgs-b + resampling
		2 threshholding
		3 direct + Metropolis sampling of the GP posetior
		4 multiple start l-bfgs-b + Metropolis sampling of the GP posetior
"""

def create_uai_tests(policy_size, nadaptations=10,
                     nruns_per_param_update=100, nsamples=10**4, strategy = 3, exp_temp = 10, MH_rate = 0.3**2):
    
    common_args = \
        [policy_size, nadaptations, nruns_per_param_update, nsamples, strategy, exp_temp, MH_rate]

    tests = \
        map(lambda x:
            (x[0], create_tests_for_ising_model(*(x[1] + common_args))),
            [
	    #('ferro2d', [ferro2d, (90, 90), (1, 300)]),
             ('isg3d', [isg3d, (1, 25), (1, 50)])#,
             #('rbm', [rbm, (1, 20), (1, 50)])
             ])
    
    return dict(tests)

default_policy_size = 1000
uai_tests = create_uai_tests(default_policy_size)

uai_kawasaki_tests = map(lambda (k, v): v['kawasaki'], uai_tests.iteritems())
uai_im_unif_tests = map(lambda (k, v): v['im_unif'], uai_tests.iteritems())
uai_im_expert_tests = map(lambda (k, v): v['im_expert'], uai_tests.iteritems())
uai_im_ego_tests = map(lambda (k, v): v['im_ego'], uai_tests.iteritems())


uai_im_ego_premixed_tests = map(lambda (k, v): v['im_ego_premixed'], uai_tests.iteritems())

def run_test(testtype, args, nruns, clear_all=False):
    test = testtype(*args)

    if clear_all:
        test.clear_all_results()

    for i in xrange(nruns):
        test.run()
        test.save_results()

def run_tests(tests, nruns, clear_all=False):
    map(lambda (testtype, args):
        run_test(testtype, args, nruns, clear_all),
        tests)

def instantiate_test(testdesc):
    testtype, args = testdesc
    test = testtype(*args)
    return test

import re

def plot_result(testtype, args, res_num):
    test = testtype(*args)
    test.load_results(res_num)
    
    test.plot_E_samples()
    if re.match('IM_EGO', test.name):
	test.plot_gp_mean()
	test.plot_policy_as_histogram()
	test.plot_acf(start = 1000)
    else:
	test.plot_acf()

def plot_tests(tests, res_num = None):
    for test in tests:
	plot_result(test[0], test[1], res_num)


run_tests(uai_im_ego_tests, 1)
plot_tests(uai_im_ego_tests)
#plot_avg_energy(uai_im_ego_tests[0])