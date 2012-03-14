from __future__ import division
from numpy import *
import ctypes as c

def make_ctypes_for_graph_repr \
    (nbrs_list_i, incident_edge_list_i, edges_i):
    # Convert arguments into ctypes-compatible types.
    num_edges_c = c.c_int32(len(edges_i))
    edges_i_c = edges_i.ctypes.data_as(c.POINTER(c.c_int32))

    nbrs_list_ptrs, nbrs_list_lens = \
        zip(*map(lambda x: (x.ctypes.data, len(x)), nbrs_list_i))

    # References to these numpy arrays must be kept while the C code
    # runs, so that it does not get garbage-collected and result in
    # segmentation faults.
    nbrs_list_ptrs = array(nbrs_list_ptrs, dtype=uint)
    nbrs_list_lens = array(nbrs_list_lens, dtype=int32)

    nbrs_list_ptrs_c = \
        nbrs_list_ptrs.ctypes.data_as(c.POINTER(c.c_size_t))
    nbrs_list_lens_c = nbrs_list_lens.ctypes.data_as(c.POINTER(c.c_int32))

    incident_edge_list_ptrs, incident_edge_list_lens = \
        zip(*map(lambda x: (x.ctypes.data, len(x)),
                 incident_edge_list_i))

    # References to these numpy arrays must be kept while the C code
    # runs, so that it does not get garbage-collected and result in
    # segmentation faults.
    incident_edge_list_ptrs = \
        array(incident_edge_list_ptrs, dtype=uint)
    incident_edge_list_lens = \
        array(incident_edge_list_lens, dtype=int32)

    incident_edge_list_ptrs_c = \
        incident_edge_list_ptrs.ctypes.data_as(c.POINTER(c.c_size_t))
    incident_edge_list_lens_c = \
        incident_edge_list_lens.ctypes.data_as(c.POINTER(c.c_int32))

    num_nodes_c = c.c_int32(len(nbrs_list_i))

    return num_nodes_c, num_edges_c,\
           edges_i_c, \
           nbrs_list_ptrs, nbrs_list_ptrs_c, \
           nbrs_list_lens, nbrs_list_lens_c, \
           incident_edge_list_ptrs, incident_edge_list_ptrs_c, \
           incident_edge_list_lens, incident_edge_list_lens_c

def make_ctypes_for_graph_weights(w):
    w_c = w.ctypes.data_as(c.POINTER(c.c_double))
    return w_c

def make_ctypes_for_node_activations(S):
    S_c = S.ctypes.data_as(c.POINTER(c.c_int32))
    return S_c
