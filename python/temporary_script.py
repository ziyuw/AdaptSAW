from __future__ import division
import tests
import ising as isg

import sys 
from numpy import *
import matplotlib.pyplot as plt
import utils
import os


if len(sys.argv) > 1:
    test_name = sys.argv[1]
else:
    test_name = 'ferrof'

rewrite = True

# Run tests based on the tests in the UAI 2010 paper.
if test_name == 'ferro':
    print "Running Ferro2D_UAI."
    test = isg.Ferro2D_UAI(60, 60)
elif test_name == 'isg3d':
    print "Running ISG3D_UAI."
    name = 'isg3d'
    if not os.path.exists(name + '.mat') or rewrite:
	test = isg.ISG3D_UAI(9, 9, 9)
elif test_name == 'rbm':
    print "Running RBM_UAI."
    name = 'RBM_UAI'
    test = isg.RBM_UAI()
elif test_name == 'ferrof':
    print "Running Ferro2D_Frustrated."
    name = 'Ferro2D_Frustrated'
    if not os.path.exists(name + '.mat') or rewrite:
	test = isg.Ferro2D_Frustrated(60, 60)
elif test_name == 'chi':
    print "Running Chimera."
    name = 'Chimera'
    if not os.path.exists(name + '.mat') or rewrite:
	test = isg.Chimera(4,4,4)

import scipy.io as sio

if not os.path.exists(name + '.mat') or rewrite:
    sio.savemat(name + '.mat', {'isg_J':test.J, 'isg_h':test.h, 'isg_edges': test.edges_i, 'incident_list': test.incident_edge_list_i, 'S_in': test.get_default_initial_state(), 'beta_true':test.beta_true })

if test_name == 'ferrof':
    import tables as t
    
    mat_contents = sio.loadmat('Ferro2D_Frustrated.mat')
    
    hdf5_fn='adaptive_saw.h5'
    h5file = t.openFile(hdf5_fn, mode = "a")
    root = h5file.root

    if not h5file.__contains__('/Ferro2DF_Param'):
	group = h5file.createGroup('/', 'Ferro2DF_Param', '', None, True)
    else:
	group = h5file.getNode('/Ferro2DF_Param')

    if not h5file.__contains__('/Ferro2DF_Param/h'):
	table1 = h5file.createArray(group, "h", mat_contents['isg_h'])
    elif rewrite:
	h5file.removeNode('/Ferro2DF_Param/h')
	table1 = h5file.createArray(group, "h", mat_contents['isg_h'])
    
    if not h5file.__contains__('/Ferro2DF_Param/J'):
	table1 = h5file.createArray(group, "J", mat_contents['isg_J'])
    elif rewrite:
	h5file.removeNode('/Ferro2DF_Param/J')
	table1 = h5file.createArray(group, "J", mat_contents['isg_J'])
    
    h5file.flush()
    h5file.close()
    
if test_name == 'chi':
    import tables as t
    
    mat_contents = sio.loadmat('Chimera.mat')
    
    hdf5_fn='adaptive_saw.h5'
    h5file = t.openFile(hdf5_fn, mode = "a")
    
    root = h5file.root

    if not h5file.__contains__('/Chimera_Param'):
	group = h5file.createGroup('/', 'Chimera_Param', '', None, True)
    else:
	group = h5file.getNode('/Chimera_Param')

    if not h5file.__contains__('/Chimera_Param/h'):
	table1 = h5file.createArray(group, "h", mat_contents['isg_h'])
    elif rewrite:
	h5file.removeNode('/Chimera_Param/h')
	table1 = h5file.createArray(group, "h", mat_contents['isg_h'])
    
    if not h5file.__contains__('/Chimera_Param/J'):
	table1 = h5file.createArray(group, "J", mat_contents['isg_J'])
    elif rewrite:
	h5file.removeNode('/Chimera_Param/J')
	table1 = h5file.createArray(group, "J", mat_contents['isg_J'])
    
    h5file.flush()
    h5file.close()