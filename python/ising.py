from __future__ import division
from numpy import *

import graphs as g
import tables as t

def get_effective_fields(edges_i, incident_edge_list_i, J, h, S):
    # This function computes h_effective, the "effective" local field
    # at each site. Includes the external field h. h_effective can be
    # used to get the energy in a cheap way:
    # E = -((h_effective-h)'*S/2 + h'*S)
    h_effective = array(h, dtype=double)
    num_nodes = h.shape[0]

    for curr_node in xrange(num_nodes):
        incident_edges = incident_edge_list_i[curr_node]
        num_nbrs = len(incident_edges)

        for curr_edge in incident_edges:
            v1 = edges_i[curr_edge, 0]
            v2 = edges_i[curr_edge, 1]

            curr_nbr = None
            if v1 == curr_node:
                # v2 is the neighbor of curr_node
                curr_nbr = v2
            else:
                # v1 is the neighbor of curr_node
                curr_nbr = v1

            # update the effective field
            h_effective[curr_node] += \
                double(S[curr_nbr]) * J[curr_edge];

    return h_effective

from graphs import make_incident_edge_list_from_edges
from graphs import make_edges_from_nbrs_list

def rbm_to_ising(W, b_visible, b_hidden):
    # b_hidden is a vector of biases for the hidden units.
    # b_visible is a vector of biases for the visible units.
    # W is the n_visible X n_hidden weight matrix.
    
    # Boltzmann machines have a {0, 1} state space and Ising models
    # have a {-1, 1} state space. The input weights are with respect
    # to the Boltzmann machine state space and the output weights are
    # with respect to the Ising model state space.

    n_visible, n_hidden = W.shape
    num_nodes = n_visible + n_hidden
    num_edges = n_visible * n_hidden # bipartite

    visible_nodes = range(n_visible)
    hidden_nodes = map(lambda x: x + n_visible, xrange(n_hidden))

    nbrs_list_i = []
    for visible_node_idx in visible_nodes:
        nbrs_list_i.append(hidden_nodes)
    for hidden_node_idx in hidden_nodes:
        nbrs_list_i.append(visible_nodes)
        
    nbrs_list_i = map(lambda x: array(x, dtype=uint).transpose(),
                      nbrs_list_i)
    
    # Construct incident edge list and list of edges
    edges_i = make_edges_from_nbrs_list(nbrs_list_i)
    incident_edge_list_i = \
        make_incident_edge_list_from_edges(edges_i, num_nodes)

    J = zeros([num_edges, 1])
    h = zeros([num_nodes, 1])

    for i in xrange(n_visible):
        h[i] = b_visible[i] / 2 + sum(W[i, :]) / 4

    for i in xrange(n_hidden):
        h[i + n_visible] = b_hidden[i] / 2 + sum(W[:, i]) / 4

    # Same order as edges were labeled. But this really only depends
    # on the fact that the edge[0]-endpoint is a visible node and
    # edge[1]-endpoint is a hidden node.
    curr_edge = 0
    for edge in edges_i:
        J[curr_edge] = W[edge[0], edge[1] - n_visible] / 4
        curr_edge += 1

    return edges_i, nbrs_list_i, incident_edge_list_i, J, h

class IsingModel(object):
    name = None
    
    beta_true = None # beta to actually use
    beta_crit = None # freezing beta
    J = None # pairwise potentials
    h = None # biases

    edges_i = None
    nbrs_list_i = None
    incident_edge_list_i = None

    saved_res_fields = set(['beta_true', 'beta_crit'])

    def get_ground_state(self):
        S_grnd = ones([self.get_num_nodes(), 1], dtype=int32)
        return S_grnd

    # Default reference state is ground state
    def get_default_reference_state(self):
        S_ref = self.get_ground_state()
        return S_ref

    # Default initial state is reference state with half of its bits flipped
    def get_default_initial_state(self):
        S_in = array(self.get_default_reference_state())
        S_in[0:self.get_num_nodes()/2] = -1
        random.shuffle(S_in)

        return S_in

    def get_num_edges(self):
        return len(self.edges_i)

    def get_num_nodes(self):
        return len(self.nbrs_list_i)

    def get_effective_fields(self, S_in):
        h_effective = \
            get_effective_fields\
            (self.edges_i, self.incident_edge_list_i, \
             self.J, self.h, S_in)
        return h_effective

    def compute_energy(self, S_in, h_effective=None):
        if h_effective is None:
            h_effective = self.get_effective_fields(S_in)
        
        E = dot(-double(S_in.transpose()), self.h) - \
            dot(double(S_in.transpose()),
                h_effective - self.h) / 2

        return E

class GridISG2D(IsingModel):
    name = "GridISG2D"
    beta_true = 1/2.27
    beta_crit = beta_true

    L = None
    M = None
    
    saved_res_fields = \
        set(['L', 'M']) | IsingModel.saved_res_fields

    def __init__(self, L, M, J, h):
        self.L = L
        self.M = M

        self.edges_i, self.nbrs_list_i, self.incident_edge_list_i = \
            g.make_periodic_grid_graph(self.L, self.M)

        self.J = J
        self.h = h

class Ferro2D_UAI(GridISG2D):
    name = "Ferro2D_UAI"
    
    def __init__(self, L, M, beta_true=None, name=None):
        super(Ferro2D_UAI, self).__init__(L, M, None, None)
        self.h = zeros([self.get_num_nodes(), 1], dtype=double)
        self.J = ones([self.get_num_edges(), 1], dtype=double)

        if beta_true is not None:
            self.beta_true = beta_true
        if name is not None:
            self.name = name
	# NOTE: This should be changed back
	#self.beta_true = 1.0/2.27

class Ferro2D_Frustrated(GridISG2D):
    name = "Ferro2D_Frustrated"

    def __init__(self, L, M, beta_true=None, name=None, hdf5_fn='adaptive_saw.h5'):
        super(Ferro2D_Frustrated, self).__init__(L, M, None, None)
        self.h = ones([self.get_num_nodes(), 1], dtype=double)
        
	h5file = t.openFile(hdf5_fn)
	try:
	    self.J = h5file.root.Ferro2DF_Param.J[:]
	    self.h = h5file.root.Ferro2DF_Param.h[:]
	except:
	    self.h = ones([self.get_num_nodes(), 1], dtype=double)
	    self.h[1:len(self.h)/2] = self.h[1:len(self.h)/2]*(-1)
	    random.shuffle(self.h)
	
	    self.J = ones([self.get_num_edges(), 1], dtype=double)
	    self.J[1:len(self.J)*0.1] = self.J[1:len(self.J)*0.1]*(-1)
	    random.shuffle(self.J)


	#NOTE: What should I set for beta_true?
        if beta_true is not None:
            self.beta_true = beta_true
        if name is not None:
            self.name = name
            
	self.beta_true = 0.8
	
	h5file.close()
	h5file = None
	

class Chimera(IsingModel):
    name = "Chimera_Graph"
    beta_true = 1.0
    beta_crit = beta_true
    
    saved_res_fields = \
        IsingModel.saved_res_fields

    def __init__(self, m, n, l, J = None, h = None, hdf5_fn='adaptive_saw.h5'):

        self.edges_i, self.nbrs_list_i, self.incident_edge_list_i, self.num_nodes = \
            g.make_chimera_graph(m,n,l)

	self.beta_true = 1.0
	self.name = "Chimera_Graph"
	
	h5file = t.openFile(hdf5_fn)
	try:
	    self.J = h5file.root.Chimera_Param.J[:]
	    self.h = h5file.root.Chimera_Param.h[:]
	except:
	    self.h = zeros([self.get_num_nodes(), 1], dtype=double)
	    from numpy.random import rand
	    self.J = rand(self.get_num_edges(), 1)*2-1
	h5file.close()
	h5file = None


class Fully_Connected_Graph(IsingModel):
    name = "Fully_Connected_Graph"
    beta_true = 1.0
    beta_crit = beta_true
    
    saved_res_fields = \
        IsingModel.saved_res_fields

    def __init__(self, num_nodes, J, h):
	self.num_nodes = num_nodes

        self.edges_i, self.nbrs_list_i, self.incident_edge_list_i = \
            g.make_fully_connected_graph(self.num_nodes)

	self.beta_true = 1.0
	self.name = "Fully_Connected_Graph"
	
	if J is None and h is None:
	    self.h = ones([self.get_num_nodes(), 1], dtype=double)
	    from numpy.random import rand
	    self.J = rand(self.get_num_edges(), 1)*2-1
	else:
	    self.J = J
	    self.h = h


class CubeISG3D(IsingModel):
    name = "CubeISG3D"
    beta_true = 1.0
    #beta_true = 0.8
    beta_crit = 0.9

    L = None
    M = None
    N = None

    saved_res_fields = \
        set(['L', 'M', 'N']) | IsingModel.saved_res_fields

    def __init__(self, L, M, N, J, h):
        self.L = L
        self.M = M
        self.N = N

        self.edges_i, self.nbrs_list_i, self.incident_edge_list_i = \
            g.make_periodic_cube_graph(L, M, N)

        self.J = J
        self.h = h

class ISG3D_UAI(CubeISG3D):
    name = "ISG3D_UAI"
    
    def __init__(self, L, M, N, hdf5_fn='adaptive_saw.h5'):
        J, h = self.load_potentials(hdf5_fn)
        super(ISG3D_UAI, self).__init__(L, M, N, J, h)

    def load_potentials(self, hdf5_fn):
        h5file = t.openFile(hdf5_fn)
        J = h5file.root.ISG3D.J[:]
        h = h5file.root.ISG3D.h[:]
        h5file.close()
        h5file = None

        return J, h

class RBM(IsingModel):
    name = "RBM"
    beta_true = 1.0
    #beta_true = 0.8

    # W - pairwise-potentials, n_visible x n_hidden
    # b - biases for hidden layer    
    # c - biases for visible layer
    def __init__(self, W, b, c):
        self.edges_i, self.nbrs_list_i, self.incident_edge_list_i, \
        self.J, self.h \
        = rbm_to_ising(W, c, b)

    # Default initial state of RBM is reference state with two-thirds
    # of its bits flipped
    def get_default_initial_state(self):
        S_in = array(self.get_default_reference_state())
        S_in[0:2*self.get_num_nodes()/3] = -1
        random.shuffle(S_in)

        return S_in

class RBM_UAI(RBM):
    name = "RBM_UAI"

    def __init__(self, hdf5_fn='adaptive_saw.h5'):
        W, b, c = self.load_potentials(hdf5_fn)
        super(RBM_UAI, self).__init__(W, b, c)

    def load_potentials(self, hdf5_fn):
        h5file = t.openFile(hdf5_fn)
        W = h5file.root.RBM.W[:].transpose()
        b = h5file.root.RBM.b[:]
        c = h5file.root.RBM.c[:]

        h5file.close()
        h5file = None

        return W, b, c

class RBM_UAI_BLOCK(RBM):
    name = "RBM_UAI_BLOCK"

    def __init__(self, hdf5_fn='adaptive_saw.h5'):
        W, b, c = self.load_potentials(hdf5_fn)
        self.W = W
        self.b = b
        self.c = c

    def load_potentials(self, hdf5_fn):
        h5file = t.openFile(hdf5_fn)
        W = h5file.root.RBM.W[:].transpose()
        b = h5file.root.RBM.b[:]
        c = h5file.root.RBM.c[:]

        h5file.close()
        h5file = None

        return W, b, c