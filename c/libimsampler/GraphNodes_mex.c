#include "GraphNodes.h"
#include "GraphNodes_mex.h"

GraphNodes* build_graphnodes_mex
    (const mxArray* mx_p_nbrs_list,
     const mxArray* mx_p_incident_edge_list,
     int num_nodes)
{
    GraphNodes *graph_nodes = (GraphNodes*)malloc(num_nodes*sizeof(GraphNodes));

    for (int i = 0; i < num_nodes; ++i)
    {
        const int *dims = mxGetDimensions(mxGetCell(mx_p_nbrs_list,i));
        int num_nbrs = dims[0];
        graph_nodes[i].num_nbrs = num_nbrs;
        graph_nodes[i].nbrs_list = (int*)mxGetData(mxGetCell(mx_p_nbrs_list,i));
        graph_nodes[i].incident_edge_list = 
            (int*)mxGetData(mxGetCell(mx_p_incident_edge_list,i));
    }

    return graph_nodes;
}
int** build_graphedges_mex(const int* p_edges, int num_edges)
{
    int** graph_edges = (int**)malloc(num_edges*sizeof(int*));

    for (int i = 0; i < num_edges; ++i)
    {
        graph_edges[i] = (int*)malloc(2*sizeof(int));
        /* Note the Fortran-Style storage */
        int v1 = p_edges[i];
        int v2 = p_edges[i+num_edges];
        graph_edges[i][0] = v1;
        graph_edges[i][1] = v2;
    }

    return graph_edges;
}

