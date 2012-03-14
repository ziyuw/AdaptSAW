#include "GraphNodes.h"
#include <stdlib.h>

void destroy_graphnodes(const GraphNodes* graph_nodes)
{
    /* Note we DO NOT free the lists of graph_nodes here as they
       "belong" to Matlab; see the assignment above.*/
    free((GraphNodes *)graph_nodes);
}
void destroy_graphedges(const int** graph_edges, int num_edges)
{
    /* FREE: */
    for (int i = 0; i < num_edges; ++i)
    {
        free((int *)graph_edges[i]);
    }
    free(graph_edges);
}

