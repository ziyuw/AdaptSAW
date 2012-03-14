#ifndef __SAWHEAPS_H_
#define __SAWHEAPS_H_

/* 15.Nov.07 */

/* Header file for the structures and functions used in the heap
   implementation of the SAW simulator. */

/* Interesting tidbit: given two properly declared and/or allocated
   structures, I have tested that if we assign one to the other, all
   the elements are copied. 

   BE CAREFUL, though: if the structs contain "arrays" i.e. ptrs, you
   may lose the pointer if you're careless!


*/

/* Individual "cells" in the heap */
typedef struct { 
  double weight; /* The key */
  double sum_curr_and_below; /* Sum of current key and all weights
				"below" this heap element.  */
  int graph_node_index; /* The graph node corresponding to
			   this heap element */
} SAWHeapElement;

typedef SAWHeapElement* pSAWHeapElement;

/* The heap */
typedef struct {
  
  pSAWHeapElement SAWHeapArray;   /* Array of SAW Heap elements... */

  int numSAWHeapElements; /* Number of entries in above array */

  int maxSAWHeapElements; /* Maximum number of elements the heap may have (num_nodes, say.) */

  int* pReverseNodeIndex; /* The array of length maxSAWHeapElements:
			     entry i says where "node" i lies in
			     SAWHeapArray. Used for updating
			     neighbors. If the node is not in, entry
			     is -1. */
} SAWHeap;

typedef SAWHeap* pSAWHeap;


pSAWHeap CreateSAWHeap( int maxSAWHeapElements );
pSAWHeap DestroySAWHeap( pSAWHeap p_saw_heap );
void AddToSAWHeap( pSAWHeap p_saw_heap, double new_item_weight, int new_item_node_index );
void removeFromSAWHeap( pSAWHeap p_saw_heap, int saw_heap_array_ind );
void changeSAWHeapElementKey( pSAWHeap p_saw_heap, int saw_heap_array_ind, double new_key );
int sampleFromSAWHeap( pSAWHeap p_saw_heap );

#endif // __SAWHEAPS_H_

