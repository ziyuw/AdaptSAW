#ifndef __GRAPHNODES_CTYPES_H__
#define __GRAPHNODES_CTYPES_H__

#include "GraphNodes.h"
#include <stdlib.h>

GraphNodes* build_graphnodes_ctypes
    (int num_nodes,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens);
int** build_graphedges_ctypes(int *edges, int num_edges);

#endif // __GRAPHNODES_CTYPES_H__

