#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* Use -I ~/project/MersenneTwister */
#include <mt19937ar.h>    
#include "SAWHeaps.h"


/* 20.Nov.07: ANNOTATE THIS STUFF! ESP THE SAMPLING AND ADD/REMOVE FUNCTIONS! */

/* Declarations of functions local to this file: */
static void updateSAWHeapWeightSums(pSAWHeap p_saw_heap, int saw_heap_array_el);
static double calculateWeightSum(pSAWHeap p_saw_heap, int elem_ind);

/* Re: the heap implementation. I have not modified the heap array to
   return something special if the heap is empty (one way is to have
   the first element always be -Inf) 
   
   Right now, check for emptiness by seeing if numSAWHeapElements ==
   0. If it is, all entries in the array can be assumed to be
   nonsense */

pSAWHeap CreateSAWHeap(int maxSAWHeapElements)
{
    /* Function to create the SAW heap. Returns a pointer to the structure.
       Inside: 1) The SAW Heap array is allocated to maxSAWHeapElements
               2) The pReverseNodeIndex is set to all -1.
     */

    int i;
    pSAWHeap p_saw_heap;

    p_saw_heap = (pSAWHeap)malloc( sizeof(SAWHeap) );

    p_saw_heap->SAWHeapArray = (pSAWHeapElement)malloc( sizeof(SAWHeapElement)*maxSAWHeapElements );

    /* Initially empty: */
    p_saw_heap->numSAWHeapElements = 0;

    p_saw_heap->maxSAWHeapElements = maxSAWHeapElements;

    p_saw_heap->pReverseNodeIndex = (int*)malloc( maxSAWHeapElements*sizeof(int) );

    for (i=0; i<maxSAWHeapElements; i++) {
        (p_saw_heap->pReverseNodeIndex)[i] = -1;
    }

    return p_saw_heap;
} 

pSAWHeap DestroySAWHeap(pSAWHeap p_saw_heap)
{
    /* Wipe out the heap structure.
       Note this has nothing to do with whether the heap is "empty" or
       not; the arrays have been allocated nontheless (of course the entries are
       undefined in that case.)
     */

    free( p_saw_heap->pReverseNodeIndex );
    free( p_saw_heap->SAWHeapArray  );
    free( p_saw_heap );

    return( (pSAWHeap)NULL );
}


void AddToSAWHeap(pSAWHeap p_saw_heap,
                  double new_item_weight, int new_item_node_index)
{
    /* Function to add a new item to the SAW Heap with weight
       new_item_weight. The index in the graph corresponding to this
       node is new_item_node_index.

       Also the pReverseIndex of the node number is set to the entry in
       the heap array containing the new element.

       After performing the addition and "heapifying" with the usual
       bottom-up propagation, the sums of "weights below" are updated to
       the root. 
     */

    int last_el_ind, curr_el_ind, parent_el_ind, parent_graph_node;

    /* Increase size of the heap, and last element is 1 minus the new size: */
    last_el_ind = ++(p_saw_heap->numSAWHeapElements) - 1;
    curr_el_ind = last_el_ind; 

    /* This jumps right to the end if the heap was initially empty */
    while ( curr_el_ind > 0 ) {

        parent_el_ind = (int)( ( curr_el_ind + 1)/2 )  - 1;

        if ( ((p_saw_heap->SAWHeapArray)[parent_el_ind]).weight < new_item_weight ) {

            /* First update the array of reverse node indices: */
            parent_graph_node = ( (p_saw_heap->SAWHeapArray)[parent_el_ind] ).graph_node_index;
            (p_saw_heap->pReverseNodeIndex)[ parent_graph_node ] = curr_el_ind;

            /* A pointwise assignment of the struct... Note that at this point
               the sum_curr_and_below field is wrong. Also this is one-way... */
            (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[parent_el_ind];

            curr_el_ind = parent_el_ind;

        } else {
            /* We're outta here! */
            break;
        }

    } /* end while */

    /* Will jump right here if the heap was initially empty. First,
       assign the fields (except the sum...) */
    ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).weight = new_item_weight;
    ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).graph_node_index = new_item_node_index;

    /* Then the reverse index*/
    (p_saw_heap->pReverseNodeIndex)[ new_item_node_index ] = curr_el_ind;

    /* Now adjust all sums along a path from the new last element to the root: */
    updateSAWHeapWeightSums( p_saw_heap, (p_saw_heap->numSAWHeapElements)-1 );

    return;
}


void removeFromSAWHeap(pSAWHeap p_saw_heap, int saw_heap_array_ind)
{
    /* Function to remove element saw_heap_array_el from the heap array.
       Special case: if it's zero, it's the "usual" delete-max
       operation. 

       NOTE: Trying to delete an element that isnt in the heap is
             undefined. 

       Elegant way it works (stole the idea from various places on the web:)

       1) Set the key of the target element to Inf, then bubble it up to the
          root.

       2) Delete the max element of the heap.
     */

    int curr_el_ind, parent_el_ind, left_child_ind, right_child_ind, P1, parent_graph_node, child_graph_node, doomed_graph_node, parent_of_last, child_ind, j;
    /* The last element prior to decreasing the size: */
    SAWHeapElement last_saw_heap_el;

    curr_el_ind = saw_heap_array_ind ;

    /* The graph node of the heap array entry that will be trashed: */
    doomed_graph_node = ( (p_saw_heap->SAWHeapArray)[saw_heap_array_ind] ).graph_node_index;

    /* while the current element has parents: */
    while ( curr_el_ind > 0 ) {

        parent_el_ind = (int)( ( curr_el_ind + 1)/2 )  - 1;

        /* No need to check dominance over parents; it's guaranteed to be
           bigger...  */

        /* First update the array of reverse node indices: */
        parent_graph_node = ( (p_saw_heap->SAWHeapArray)[parent_el_ind] ).graph_node_index;
        (p_saw_heap->pReverseNodeIndex)[ parent_graph_node ] = curr_el_ind;

        /* A pointwise assignment of the struct... Note that at this point
           the sum_curr_and_below field is wrong. Will be fixed... */
        (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[parent_el_ind];
        /* Fix the sum field: */
        ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).sum_curr_and_below = calculateWeightSum( p_saw_heap, curr_el_ind );

        curr_el_ind = parent_el_ind;
    }

    /* Set the key of the current element (should be root) to largest
       possible. All other fields are undefined; the element is about to be
       trashed anyway. */
    ((p_saw_heap->SAWHeapArray)[ curr_el_ind ]).weight = __builtin_inf();

    /* Now we delete the largest element (root.) */

    /* First, store the last element (prior to reducing the size) */
    last_saw_heap_el = (p_saw_heap->SAWHeapArray)[ p_saw_heap->numSAWHeapElements - 1 ];

    /* P1 is the parent of the (prior to reducing size) last node. We
       store it to later adjust the weight sums from it to the root. */
    P1 = (int)( p_saw_heap->numSAWHeapElements / 2 )  - 1;

    /* Now reduce the size: */
    (p_saw_heap->numSAWHeapElements)--;

    /* Start at root: */
    curr_el_ind = 0;

    /* Loop until at most the node which is the parent of the last node
       in the reduced system. Note that this ensures that upon entry to
       the while, we will always have at least the left child. */
    parent_of_last = (int)( p_saw_heap->numSAWHeapElements / 2 )  - 1;
    while ( curr_el_ind <= parent_of_last ) {

        left_child_ind = 2*curr_el_ind + 1;
        right_child_ind = left_child_ind  + 1; 


        /* 2 cases handled separately: Either 2 children or 1 (the left) child.
         */
        if ( right_child_ind < p_saw_heap->numSAWHeapElements ) {
            /* 2 children case. To ensure that it stays a heap, the path
               needs to go down the larger of the 2 children, if one is
               larger. Otherwise, choose the left one arbitrarily. */

            if ( ((p_saw_heap->SAWHeapArray)[ left_child_ind ]).weight >= ((p_saw_heap->SAWHeapArray)[ right_child_ind ]).weight ) {
                child_ind = left_child_ind;
            } else {
                child_ind = right_child_ind;
            }

            /* Now comparisons are done with the "working" child... */
            if ( ((p_saw_heap->SAWHeapArray)[ child_ind ]).weight > last_saw_heap_el.weight ) {
                /* curr el gets the info of left child */

                /* First update the array of reverse node indices: */
                child_graph_node = ( (p_saw_heap->SAWHeapArray)[child_ind] ).graph_node_index;
                (p_saw_heap->pReverseNodeIndex)[ child_graph_node  ] = curr_el_ind;

                /* A pointwise assignment of the struct... Note that at this point
                   the sum_curr_and_below field is wrong. */
                (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[child_ind];

                curr_el_ind = child_ind; 

            } else { /* a valid heap; get outta here */
                break;
            }

        } else {
            /* Only have left child. Must be a left one, i.e. no need to
               figure out which is the bigger one */

            if ( ((p_saw_heap->SAWHeapArray)[ left_child_ind ]).weight > last_saw_heap_el.weight ) {
                /* curr el gets the info of left child */

                /* First update the array of reverse node indices: */
                child_graph_node = ( (p_saw_heap->SAWHeapArray)[left_child_ind] ).graph_node_index;
                (p_saw_heap->pReverseNodeIndex)[ child_graph_node ] = curr_el_ind;

                /* A pointwise assignment of the struct... Note that at this point
                   the sum_curr_and_below field is wrong. */
                (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[left_child_ind];

                curr_el_ind = left_child_ind;

            } else { /* a valid heap; get outta here */
                break;
            }
        }
    } /* end while */

    /* Finally, assign curr_el's info to the stuff in the last_saw_heap_el */
    /* Then, fix the sums: 
       1) from P1 to root.
       2) from curr_el_ind to root.
     */
    /* Only do this if the heap is nonzero size now. IS THIS CONDITION
       needed?!  */
    if ( p_saw_heap->numSAWHeapElements > 0 ) {
        ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).weight = last_saw_heap_el.weight;
        ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).graph_node_index = last_saw_heap_el.graph_node_index;
        /* Don't forget to set the index graph index of the current
           node!!!!  */
        (p_saw_heap->pReverseNodeIndex)[last_saw_heap_el.graph_node_index] = curr_el_ind;

    }

    /* The reverse index of the trashed node is simply set to -1 */
    (p_saw_heap->pReverseNodeIndex)[ doomed_graph_node ] = -1;

    /* First update along P1 to root */
    updateSAWHeapWeightSums( p_saw_heap, P1 );
    /* Then, along curr_el_ind to root */
    updateSAWHeapWeightSums( p_saw_heap, curr_el_ind );

    return;
}

void changeSAWHeapElementKey
    (pSAWHeap p_saw_heap, int saw_heap_array_ind, double new_key)
{
    /* 21.Nov.07: Function to change the key of saw_heap_array_ind to
       new_key (keeping all other data the same) and keeping it a valid
       heap. 

       Behavior undefined if you try to change the key of an element
       past the current heap size!

       Current way it works: we put two loops after each other, one to
       propagate to the root (i.e. new_key of the element is larger than
       its parent's key ) and one to propagate down to the children.

       Rationale: if it is larger than its parent, it will become a
       valid heap after leaving the while. It enters the while checking
       if it's smaller than the children, but will break leave right
       away without changing anything.

       In contrast, if it's smaller than one of its children, then
       (since it was a valid heap before the new key) it cannot be
       larger than the parent. The parent loop thus breaks with no
       modification to anything, and we pass to the child loop. This
       loop then makes it a proper heap by passing it down.

       21.Nov.07: Appears to work fine.
     */

    int curr_el_ind = saw_heap_array_ind, parent_el_ind, left_child_ind, right_child_ind, child_ind, heap_size, graph_node, child_graph_node, parent_graph_node;

    heap_size = p_saw_heap->numSAWHeapElements;

    /* The graph node of the re-keyed element will not be changed; only the key is. So store it
       here till the end: */
    graph_node = ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).graph_node_index;

    /* First loop: only enter if curr_el_ind has a parent */
    while ( curr_el_ind > 0) {

        parent_el_ind = (int)( ( curr_el_ind + 1)/2 )  - 1;

        /* New key larger than parent of current element's? */
        if ( ((p_saw_heap->SAWHeapArray)[parent_el_ind]).weight < new_key ) {
            /* Yes... */

            /* First update the array of reverse node indices: */
            parent_graph_node = ( (p_saw_heap->SAWHeapArray)[parent_el_ind] ).graph_node_index;
            (p_saw_heap->pReverseNodeIndex)[ parent_graph_node ] = curr_el_ind;

            /* A pointwise assignment of the struct... Note that at this point
               the sum_curr_and_below field is wrong. Will be fixed... */
            (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[parent_el_ind];
            /* Fix the sum field: */
            ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).sum_curr_and_below = calculateWeightSum( p_saw_heap, curr_el_ind );

            curr_el_ind = parent_el_ind;

        } else {
            /* Note this does NOT necessarily mean we yet have a valid heap! */
            break;
        }

    } /* end while */

    /* Second loop: only enter if curr_el_ind has at least a left child... */
    /* Note this is an equivalent condition to the one in the remove
       function- work it out with the floors and <= vs < */
    while ( (2*curr_el_ind + 1) < heap_size ) {

        left_child_ind = 2*curr_el_ind + 1;
        right_child_ind = left_child_ind  + 1; 

        /* 2 cases handled separately: Either 2 children or 1 (the left) child.
         */
        if ( right_child_ind < heap_size ) {
            /* 2 children case. To ensure that it stays a heap, the path
               needs to go down the larger of the 2 children, if one is
               larger. Otherwise, choose the left one arbitrarily. */

            if ( ((p_saw_heap->SAWHeapArray)[ left_child_ind ]).weight >= ((p_saw_heap->SAWHeapArray)[ right_child_ind ]).weight ) {
                child_ind = left_child_ind;
            } else {
                child_ind = right_child_ind;
            }

            /* Now comparisons are done with the "working" child... */
            if ( ((p_saw_heap->SAWHeapArray)[ child_ind ]).weight > new_key ) {
                /* curr el gets the info of left child */

                /* First update the array of reverse node indices: */
                child_graph_node = ( (p_saw_heap->SAWHeapArray)[child_ind] ).graph_node_index;
                (p_saw_heap->pReverseNodeIndex)[ child_graph_node  ] = curr_el_ind;

                /* A pointwise assignment of the struct... Note that at this point
                   the sum_curr_and_below field is wrong. */
                (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[child_ind];

                curr_el_ind = child_ind; 

            } else { /* by now, a valid heap; get outta here */
                break;
            }

        } else {
            /* Only have left child. Must be a left one, i.e. no need to
               figure out which is the bigger one */

            if ( ((p_saw_heap->SAWHeapArray)[ left_child_ind ]).weight > new_key ) {
                /* curr el gets the info of left child */

                /* First update the array of reverse node indices: */
                child_graph_node = ( (p_saw_heap->SAWHeapArray)[left_child_ind] ).graph_node_index;
                (p_saw_heap->pReverseNodeIndex)[ child_graph_node ] = curr_el_ind;

                /* A pointwise assignment of the struct... Note that at this point
                   the sum_curr_and_below field is wrong. */
                (p_saw_heap->SAWHeapArray)[curr_el_ind] = (p_saw_heap->SAWHeapArray)[left_child_ind];

                curr_el_ind = left_child_ind;

            } else { /* by now a valid heap; get outta here */
                break;
            }
        }
    } /* end while */

    /* Finally, assign curr_el's info to the stuff in the
       last_saw_heap_el

       Then, fix the sum curr_el_ind to root.
     */
    /* Only do this if the heap is nonzero size now. IS THIS CONDITION
       needed?!  */
    if ( heap_size > 0 ) {
        ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).weight = new_key;
        ( (p_saw_heap->SAWHeapArray)[curr_el_ind] ).graph_node_index = graph_node;
        /* Don't forget to set the index graph index of the current
           node!!!!  */
        (p_saw_heap->pReverseNodeIndex)[graph_node] = curr_el_ind;

    }

    /* Update from curr_el_ind to the root */
    updateSAWHeapWeightSums( p_saw_heap, curr_el_ind );

    return;
}

int sampleFromSAWHeap(pSAWHeap p_saw_heap)
{
    /* Function to sample an element from the SAW heap based on the
       weights at the elements. Undefined behaviour if you try on an empty
       heap.

       Returned: the position in the SAWHeapArray of the sampled element. NOT
       the graph node! You can recover all the necessary info such as the node
       in the graph it corresponded to by accessing the fields of that entry of
       the array.

       Works by refining the interval that random variable R could have been
       in. Pretty straightforward.

       Note: I could have made it conceptually more simple by repeatedly
       drawing random samples to choose a path along the heap, but I did it
       this way for the challenge. Only one random sample is drawn per call to
       this function.
     */
    int curr_element = 0, num_elements, left_child, right_child;
    double interval_start = 0.0, R, total_weight, W_e, S_L;

    num_elements = p_saw_heap->numSAWHeapElements;

    /* Sample R ~ U[0,1] */
    R = genrand_real2();

    /* Scale it to lie in U[0, total_weight] */
    total_weight = ( (p_saw_heap->SAWHeapArray)[0] ).sum_curr_and_below;
    R *= total_weight;

    /* Loop while curr_element has children (i.e. at least a left child) */
    while ( (2*curr_element+1) < num_elements ) {

        /* No need to check for number of children. See comment later. */

        left_child = 2*curr_element+1;
        /* may not be "legitimate" but so what. */
        right_child = left_child + 1;

        /* weight of current el: */
        W_e = ( (p_saw_heap->SAWHeapArray)[curr_element] ).weight;

        /* sum of weights below the left child: */
        S_L = ( (p_saw_heap->SAWHeapArray)[left_child] ).sum_curr_and_below;
        /* We have no need for the corresponding right child weight sum */

        /* If the weight of the curr_el has won:
           R in [interval_start, W_e+interval_start] */
        /* Reminder: R has been scaled!! */
        if ((R > interval_start) && (R <= interval_start+W_e))
        {
            break;
        }
        else if ((R > interval_start+W_e) && (R <= interval_start+W_e+S_L))
        { /* weight of the left child has won: */

            interval_start += W_e;
            curr_element = left_child;
        }
        else
        { /* weight of the right child has won: */
            /* NOTE We can NEVER ENTER here if curr_element has no right
               child! Reason is that in that case the stored sum of
               curr_element gets no contribution of weight from that side!
               Thus, no need to check... */

            interval_start += ( W_e + S_L );
            curr_element = right_child;
        }

    } /* end while */

    return(curr_element);
}

static void updateSAWHeapWeightSums(pSAWHeap p_saw_heap, int saw_heap_array_el)
{
    /* Function to update the sums of weights in the SAW heap array
       from node saw_heap_array_el to the root.
     */

    int curr_el_ind = saw_heap_array_el;

    /* while current element is "valid" */
    while ( curr_el_ind >= 0 ) {

        ( (p_saw_heap->SAWHeapArray)[ curr_el_ind ] ).sum_curr_and_below = calculateWeightSum( p_saw_heap, curr_el_ind );

        /* current node becomes its parent: */
        curr_el_ind = (int)( ( curr_el_ind + 1)/2 )  - 1;

    } /* end while */
}

static double calculateWeightSum(pSAWHeap p_saw_heap, int elem_ind)
{
    /* Little bookeeping function to add the weight of element elem_ind
       in the heap to that of its two children.

       The only thing stopping this from being totally trivial is that a
       node may not actually have 2 children (in which case, the missing
       child contributes zero) so I have written this for cleanliness.
     */
    int left_child_ind, right_child_ind;
    double left_weight_sum, right_weight_sum, s;

    left_child_ind = 2*elem_ind + 1;
    right_child_ind = left_child_ind + 1;

    /* The "defaults," unless current element has children */
    left_weight_sum = 0.0;
    right_weight_sum = 0.0;
    if ( left_child_ind < p_saw_heap->numSAWHeapElements ) {

        left_weight_sum = ( (p_saw_heap->SAWHeapArray)[ left_child_ind ] ).sum_curr_and_below;

        if ( right_child_ind < p_saw_heap->numSAWHeapElements ) {
            right_weight_sum = ( (p_saw_heap->SAWHeapArray)[ right_child_ind ] ).sum_curr_and_below;
        }

    }
    s = ( (p_saw_heap->SAWHeapArray)[ elem_ind ] ).weight + left_weight_sum + right_weight_sum;

    return(s);
}
