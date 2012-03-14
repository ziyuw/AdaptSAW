from __future__ import division
import tests
import ising as isg

import numpy.random as rn
from numpy import *
import matplotlib.pyplot as plt
import utils

def sigmoid(x):
    return 1/(1+exp(-x))

def save_current_plot(graph_name , plot_name, method_name, res_num):
    import utils
    import os.path
    import matplotlib.pyplot as plt
    
    plot_dir = "plots/%s" % (graph_name)
    plot_fn = "%s/%s_%s_r%d.pdf" % (plot_dir, plot_name, method_name, res_num)

    utils.ensure_dir_exists(os.path.dirname(plot_fn))
    plt.savefig(plot_fn, format='pdf')
    print "Saved plot to %s" % (plot_fn)

class Block_RBM_Sampler:
    def __init__(self, rbm, v = None, h = None):
	self.rbm = rbm
	self.v = v
	self.h = h
	self.name = 'bGibbs'
	
	if v is None:
	    self.v = zeros_like(self.rbm.c)
	if h is None:
	    self.h = zeros_like(self.rbm.b)

    def energy(self):
	return -dot(self.v, dot(self.rbm.W, self.h)) - dot(self.rbm.b, self.h) - dot(self.rbm.c, self.v) \
	    + dot(ones(784), dot(self.rbm.W, ones(500)))*0.25 + 0.5*dot(self.rbm.c, ones(784)) + 0.5*dot(self.rbm.b, ones(500))
    
    def assign(self, zipped, h=True):
	for i in xrange(len(zipped)):
	    if h:
		if zipped[i][1] < zipped[i][0]:
		    self.h[i] = 1.0
		else:
		    self.h[i] = 0.0
	    else:
		if zipped[i][1] < zipped[i][0]:
		    self.v[i] = 1.0
		else:
		    self.v[i] = 0.0
    
    def block_sample(self):
	zipped_h = zip(sigmoid(self.rbm.beta_true*(dot(self.rbm.W.T, self.v) + self.rbm.b)), rn.rand(self.rbm.b.shape[0]))
	self.assign(zipped_h)
	zipped_v = zip(sigmoid(self.rbm.beta_true*(dot(self.rbm.W, self.h) + self.rbm.c)), rn.rand(self.rbm.b.shape[0]))
	self.assign(zipped_v, h=False)

    def iterate(self, n):
	self.n = n
	self.energies = []
	self.energies.append(self.energy())
	for i in xrange(n):
	    self.block_sample()
	    self.energies.append(self.energy())
	
	return self.energies
	
    def plots(self, max_lag = 500, res_num = 0, energies = None, save_plots = False):
	if energies is None:
	    energies = self.energies
	plt.figure()
	plt.plot(energies)
	if save_plots:
	    save_current_plot(self.rbm.name, 'E_samples', self.name+'_'+str(self.n), res_num)
	
	plt.figure()
	plt.plot(map(lambda i: utils.compute_acf(energies[1:], i), xrange(max_lag+1)))
	plt.ylim(ymax=1.0)
	plt.ylim(ymin=0.0)
	
	plt.grid(True)
	if save_plots:
	    save_current_plot(self.rbm.name, 'Acf', self.name+'_s'+str(self.n), res_num)
	
	if not save_plots:
	    plt.show()
	
    def acf_plot_multi_energies(self, energies, max_lag = 500, res_num = 0):
	total_acf = zeros(max_lag+1)
	for energy in energies:
	    total_acf += map(lambda i: utils.compute_acf(energy[1:], i), xrange(max_lag+1))
	
	total_acf = total_acf/float(len(energies))
	
	plt.figure()
	plt.plot(total_acf)
	plt.ylim(ymax=1.0)
	plt.ylim(ymin=0.0)
	
	plt.grid(True)
	save_current_plot(self.rbm.name, 'Acf_total', self.name+'_s'+str(self.n), res_num)

if  __name__ == "__main__":
    import ASAW_config, sys
    batch = False
    num_iter = 100000
    
    if len(sys.argv) > 2:
	batch = bool(sys.argv[2])
	num_iter = int(sys.argv[1])
    elif len(sys.argv) > 1:
	num_iter = int(sys.argv[1])
    
    if batch:
	rbm = isg.RBM_UAI_BLOCK()
	v = zeros(784)
	
	total_energy = []
	
	for i in xrange(5):
	    v[0:428] = ones(428)
	    bbs = Block_RBM_Sampler(rbm, v, None)
	    total_energy.append(bbs.iterate(num_iter))

	cur_res_num = ASAW_config.get_and_inc_cur_res_num('path_config.cfg')
	bbs.acf_plot_multi_energies(total_energy, res_num = cur_res_num,)
	f = open('Block_RBM_' + str(cur_res_num) + '.npy', 'w')
	save(f, total_energy)
	
	f.close()
    else:    
	rbm = isg.RBM_UAI_BLOCK()
	v = zeros(784)
	total_energy = []
	
	v[0:428] = ones(428)
	bbs = Block_RBM_Sampler(rbm, v, None)
	import time
    
	time_1 = time.time()
	bbs.iterate(num_iter)
	time_2 = time.time()
	print "Sampling lasted: ", time_2 - time_1
	
	cur_res_num = ASAW_config.get_and_inc_cur_res_num('path_config.cfg')
	bbs.plots(res_num = cur_res_num, save_plots = True)
    