#ifndef __GRAPHNODS_MEX_H__
#define __GRAPHNODS_MEX_H__

#include "GraphNodes.h"
#include "mex.h"

GraphNodes* build_graphnodes_mex
    (const mxArray* mx_p_nbrs_list,
     const mxArray* mx_p_incident_edge_list,
     int num_nodes);
int** build_graphedges_mex(const int* p_edges, int num_edges);

#endif // __GRAPHNODS_MEX_H__

