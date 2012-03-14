from __future__ import division
from numpy import *

def make_periodic_cube_graph(vert_size, horz_size, depth):
    layer_size = vert_size * horz_size
    num_nodes = layer_size * depth

    # Construct node neighbours list
    nbrs_list_i = []
    for curr_node in xrange(num_nodes):
        nbrs_list = []
        
        l = curr_node % layer_size # move into the L-M plane
        l = floor(l / horz_size)

        m = curr_node % layer_size # move into the L-M plane
        m = m % horz_size

        n = floor(curr_node / layer_size)

        # N prev
        k = (n - 1) % depth
        nbr = k * layer_size + l * horz_size + m
        nbrs_list.append(nbr)

        # L prev
        k = (l - 1) % vert_size
        nbr = n * layer_size + k * horz_size + m
        nbrs_list.append(nbr)

        # M prev
        k = (m - 1) % horz_size
        nbr = n * layer_size + l * horz_size + k
        nbrs_list.append(nbr)

        # M next
        # don't duplicate the loop-around edge for horz_size == 2
        if (m - 1) % horz_size != (m + 1) % horz_size:
            k = (m + 1) % horz_size
            nbr = n * layer_size + l * horz_size + k
            nbrs_list.append(nbr)

        # L next
        # don't duplicate the loop-around edge for vert_size == 2
        if (l - 1) % vert_size != (l + 1) % vert_size:
            k = (l + 1) % vert_size
            nbr = n * layer_size + k * horz_size + m
            nbrs_list.append(nbr)
        
        # N next
        # don't duplicate the loop-around edge for depth == 2
        if (n - 1) % depth != (n + 1) % depth:
            k = (n + 1) % depth
            nbr = k * layer_size + l * horz_size + m
            nbrs_list.append(nbr)

        nbrs_list_i.append(map(int, nbrs_list))

    nbrs_list_i = map(lambda x: array(x, dtype=int32).transpose(),
                      nbrs_list_i)    

    # Construct incident edge list and list of edges
    edges_i = make_edges_from_nbrs_list(nbrs_list_i)
    incident_edge_list_i = \
        make_incident_edge_list_from_edges(edges_i, num_nodes)

    return edges_i, nbrs_list_i, incident_edge_list_i

def make_periodic_grid_graph(vert_size, horz_size):
    num_nodes = vert_size * horz_size

    # Construct node neighbours list
    nbrs_list_i = []
    for curr_node in xrange(num_nodes):
        nbrs_list = []
        
        l = curr_node % vert_size
        m = floor(curr_node / vert_size)

        # L prev
        k = (l - 1) % vert_size
        nbr = m * vert_size + k
        nbrs_list.append(nbr)

        # M prev
        k = (m - 1) % horz_size
        nbr = k * vert_size + l
        nbrs_list.append(nbr)

        # M next
        # don't duplicate the loop-around edge for horz_size == 2
        if (m - 1) % horz_size != (m + 1) % horz_size:
            k = (m + 1) % horz_size
            nbr = k * vert_size + l
            nbrs_list.append(nbr)

        # L next
        # don't duplicate the loop-around edge for vert_size == 2
        if (l - 1) % vert_size != (l + 1) % vert_size:
            k = (l + 1) % vert_size
            nbr = m * vert_size + k
            nbrs_list.append(nbr)

        nbrs_list_i.append(map(int, nbrs_list))
        
    nbrs_list_i = map(lambda x: array(x, dtype=int32).transpose(),
                      nbrs_list_i)

    # Construct incident edge list and list of edges
    edges_i = make_edges_from_nbrs_list(nbrs_list_i)
    incident_edge_list_i = \
        make_incident_edge_list_from_edges(edges_i, num_nodes)

    return edges_i, nbrs_list_i, incident_edge_list_i

def make_chimera_graph(m, n, l):
    """
	m: number of rows
	n: number of columns
	l: full bipartite size
    """
    
    num_nodes = m*n*l*2
    
    def node_num(node):
	return (2*l*((node[0]*n)+node[1])) + node[2]*l + node[3]
    
    def rule1(n1,n2):
	if not n1[0] == n2[0]:
	    return False
	if not n1[1] == n2[1]:
	    return False
	if not abs(n1[2] - n2[2]) == 1:
	    return False
	return True
	
	
    def rule2(n1, n2):
	if not n1[2] == n2[2] or n2[2] == 1:
	    return False
	if not n1[3] == n2[3] or not n1[1] == n2[1]:
	    return False
	if not abs(n1[0] - n2[0]) == 1 :
	    return False
	return True

    def rule3(n1, n2):
	if not n1[2] == n2[2] or n2[2] == 0:
	    return False
	if not n1[3] == n2[3] or not n1[0] == n2[0]:
	    return False
	if not abs(n1[1] - n2[1]) == 1 :
	    return False
	return True
	
    nodes = []
    for i in range(m):
	for j in range(n):
	    for u in range(2):
		for k in range(l):
		    nodes.append((i,j,u,k))
    
    i = 0
    nbrs_list_i = []
    for n1 in nodes:
	nbrs_list = []
	assert(node_num(n1) == i)
	for n2 in nodes:
	    if rule1(n1,n2):
		nbrs_list.append(node_num(n2))
	    if rule2(n1, n2):
		nbrs_list.append(node_num(n2))
	    if rule3(n1, n2):
		nbrs_list.append(node_num(n2))
	i += 1
	nbrs_list_i.append(map(int, nbrs_list))
    
    nbrs_list_i = map(lambda x: array(x, dtype=int32).transpose(), nbrs_list_i)
    edges_i = make_edges_from_nbrs_list(nbrs_list_i)
    incident_edge_list_i = make_incident_edge_list_from_edges(edges_i, num_nodes)

    return edges_i, nbrs_list_i, incident_edge_list_i, num_nodes

def make_fully_connected_graph(num_nodes):
    nbrs_list_i = []
    
    for curr_node in xrange(num_nodes):
	nbrs_list = []

	for node in xrange(num_nodes):
	    if not (node == curr_node):
		nbrs_list.append(node)

	nbrs_list_i.append(map(int, nbrs_list))
    
    nbrs_list_i = map(lambda x: array(x, dtype=int32).transpose(),
                      nbrs_list_i)

    # Construct incident edge list and list of edges
    edges_i = make_edges_from_nbrs_list(nbrs_list_i)
    incident_edge_list_i = \
        make_incident_edge_list_from_edges(edges_i, num_nodes)

    return edges_i, nbrs_list_i, incident_edge_list_i

def make_edges_from_nbrs_list(nbrs_list_i):
    num_nodes = len(nbrs_list_i)
    edges_i = []
    
    for (curr_node, nbrs) in zip(xrange(num_nodes), nbrs_list_i):
        for nbr in nbrs:
            if curr_node > nbr:
                # Avoid duplicate edges
                continue

            edges_i.append([curr_node, nbr])

    edges_i = array(edges_i, dtype=int32)
    return edges_i

def make_incident_edge_list_from_edges(edges_i, num_nodes):
    incident_edge_list_i = [[] for i in xrange(num_nodes)]
    
    edge_idx = 0    
    for edge in edges_i:
        node = edge[0]
        othernode = edge[1]

        incident_edge_list_i[node].append(edge_idx)
        incident_edge_list_i[othernode].append(edge_idx)

        edge_idx += 1

    incident_edge_list_i = \
        map(lambda x: array(x, dtype=int32).transpose(),
            incident_edge_list_i)
    
    return incident_edge_list_i


if  __name__ == "__main__":
    edges_i, nbrs_list_i, incident_edge_list_i, num_nodes = make_chimera_graph(8,8,4)
    print num_nodes
    