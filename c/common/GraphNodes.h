#ifndef __GRAPHNODES_H_
#define __GRAPHNODES_H_

/* Standard typedef for a graph node: */
typedef struct gnodes {
  int num_nbrs;
  int* nbrs_list;
  int* incident_edge_list;
} GraphNodes;

void destroy_graphnodes(const GraphNodes* graph_nodes);
void destroy_graphedges(const int** graph_edges, int num_edges);

#endif // __GRAPHNODES_H_

