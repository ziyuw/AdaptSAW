#include "GraphNodes.h"
#include "GraphNodes_ctypes.h"

GraphNodes* build_graphnodes_ctypes
    (int num_nodes,
     int **nbrs_lists, int *nbrs_list_lens,
     int **incident_edge_lists, int *incident_edge_list_lens)
{
    GraphNodes *graph_nodes = (GraphNodes*)malloc(num_nodes*sizeof(GraphNodes));

    for (int i = 0; i < (int)num_nodes; ++i)
    {
        graph_nodes[i].nbrs_list = (int*)nbrs_lists[i];
        graph_nodes[i].num_nbrs = (int)nbrs_list_lens[i];

        // The number of incident edges is equal to the number of neighbours.
        graph_nodes[i].incident_edge_list = (int*)incident_edge_lists[i];
    }

   return graph_nodes;
}

int** build_graphedges_ctypes(int *edges, int num_edges)
{
    int** graph_edges = (int**)malloc(num_edges*sizeof(int*));

    for (int i = 0; i < (int)num_edges; ++i)
    {
        graph_edges[i] = (int*)malloc(2*sizeof(int));
        // C-style storage
        int v1 = (int)edges[i*2];
        int v2 = (int)edges[i*2+1];
        graph_edges[i][0] = v1;
        graph_edges[i][1] = v2;
    }

    return graph_edges;
}

