import tests
import ising as isg

from numpy import *
import numpy.random as rdm
import matplotlib.pyplot as plt
import random
import utils

class Ising_Gibbs:
    
    def __init__(self, graph, rand_states=False):
	self.dim = graph.get_num_nodes()
	self.graph = graph
	self.nbrs = self.graph.nbrs_list_i
	self.edges = self.graph.edges_i
	
	if rand_states:
	    self.states = [random.randint(0, 1)*2-1 for i in xrange(self.dim)]
	else:
	    self.states = self.graph.get_default_initial_state()
	self.J = self.graph.J
	self.h = self.graph.h
	self.compute_energy()
	self.beta = self.graph.beta_true
	self.incident_edge_list_i = self.graph.incident_edge_list_i

    def compute_energy(self):
	
	self.energy = 0;
	for i in xrange(len(self.edges)):
	    self.energy = self.energy - self.states[self.edges[i][0]]*self.states[self.edges[i][1]]*self.J[i]
	for i in xrange(self.dim):
	    self.energy = self.energy - self.states[i]*self.h[i]

    def gibbs_step(self, i):
	curr_state = self.states[i]
	proposed_state = -1*curr_state
	
	proposed_energy = self.energy + self.h[i] * (- proposed_state + curr_state)
	for j in xrange(len(self.incident_edge_list_i[i])):
	    
	    if self.edges[self.incident_edge_list_i[i][j]][0] == i:
		other_node_index = self.edges[self.incident_edge_list_i[i][j]][1]
	    else:
		other_node_index = self.edges[self.incident_edge_list_i[i][j]][0]
	    
	    proposed_energy += self.J[self.incident_edge_list_i[i][j]] * (-proposed_state + curr_state) * self.states[other_node_index]

	if random.random() < 1/(1 + exp(-self.beta*(self.energy - proposed_energy))):
	    self.states[i] = proposed_state
	    self.energy = proposed_energy
	else:
	    self.states[i] = curr_state

    def gibbs_sweep(self, num_iter):
	
	self.energies = zeros((num_iter*self.dim**2, 1))
	
	for i in xrange(num_iter):
	    print "Iteration: ", i
	    for j in xrange(self.dim**2):
		self.gibbs_step(j)
		self.energies[i*self.dim**2+j] = self.energy
    
    def gibbs_random(self, num_iter):
	self.energies = zeros((num_iter, 1))
	indices = rdm.randint(0, self.dim, num_iter);
	
	for i in xrange(num_iter):
	    self.gibbs_step(indices[i])
	    self.energies[i] = self.energy

    def compute_acf(self, max_lag = 500, subsample = 1):
	self.acf = map(lambda i: utils.compute_acf(self.energies[::subsample], i), xrange(max_lag+1))

    def plots(self):
	plt.figure()
	line = plt.plot(self.energies, '-*')
	plt.grid(True)
	#plt.figure()
	#plt.imshow(reshape(self.states, (sqrt(self.dim), sqrt(self.dim))))

	plt.figure()
	plt.plot(self.acf)
	plt.ylim(ymax=1.0)
	plt.ylim(ymin=0.0)
	plt.show()

if  __name__ == "__main__":
    import sys 
    import ASAW_config
    
    num_iter = 10000
    if len(sys.argv) > 2:
	test_name = sys.argv[1]
	num_iter = int(sys.argv[2])
    elif len(sys.argv) > 1:
	test_name = sys.argv[1]
    else:
	test_name = 'ferro'

    if test_name == 'ferro':
	print 'Running Ferro2D_UAI for {0} Iterations.'.format(num_iter)
	test = isg.Ferro2D_UAI(60, 60)
    elif test_name == 'isg3d':
	print 'Running ISG3D_UAI for {0} Iterations.'.format(num_iter)
	test = isg.ISG3D_UAI(9, 9, 9)
    elif test_name == 'rbm':
	print 'Running RBM_UAI for {0} Iterations.'.format(num_iter)
	test = isg.RBM_UAI()
    elif test_name == 'ferrof':
	print 'Running Ferro2D_Frustrated for {0} Iterations.'.format(num_iter)
	test = isg.Ferro2D_Frustrated(60, 60)
    elif test_name == 'chi':
	print "Running Chimera."
	test = isg.test = isg.Chimera(4,4,4)
    
    temp = Ising_Gibbs(test)
    
    import time
    
    time_1 = time.time()
    temp.gibbs_random(num_iter)
    time_2 = time.time()
    print "Sampling lasted: ", time_2 - time_1
    
    print 'Saving results in ' + './Regular_Gibbs_' + test_name + '_' + str(num_iter) + '.npy'
    
    cur_res_num = ASAW_config.get_and_inc_cur_res_num('path_config.cfg')
    f = open('Regular_Gibbs_' + test_name + '_' + str(num_iter) + '.npy', 'w')
    save(f, temp.energies[1:])
    f.close()
    
    temp.compute_acf(max_lag = min(num_iter-1, 500), subsample = max(num_iter/100000, 1))
    
    temp.plots()
