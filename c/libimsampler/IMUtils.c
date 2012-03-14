#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<mt19937ar.h>
#include "GraphNodes.h"
#include "IMOutData.h"
#include "SAWHeaps.h"
#include "IM.h"
#include "IMUtils.h"

static double getPFlipGivenChosen
    (const double* h_effective, const int* S, double beta, int I);
static void updateHeapAfterFlip
    (pSAWHeap p_candidates_heap, const int* S, const double* h_effective,
     const GraphNodes* graph_nodes, const int** graph_edges, int curr_node,
     double beta);
static void singleFlipHTranslate
    (const int* S, const GraphNodes* graph_nodes, const int** graph_edges,
     double* h_effective, const double* J, int num_nodes, int curr_node);

void initIMState(IMState* p_im_state, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params)
{
    // Also does the allocations
    int i;

    p_im_state->S = (int*)malloc( p_in_model->num_nodes*sizeof(int) );
    p_im_state->h_effective = (double*)malloc( p_in_model->num_nodes*sizeof(double) );

    for ( i=0; i<p_in_model->num_nodes; i++ ) { 
        p_im_state->S[i] = p_sim_params->S_in[i];
        p_im_state->h_effective[i] = p_sim_params->h_effective_in[i];
    }

    p_im_state->E = p_sim_params->E_in;

    return;
}

void setIMState(IMState* p_im_state,
                const IMInputModel* p_in_model,
                const int *S_new, const double* h_effective_new)
{
    printf("c\n");
    for (int i = 0; i < p_in_model->num_nodes; i++)
    { 
        p_im_state->S[i] = S_new[i];
        p_im_state->h_effective[i] = h_effective_new[i];
    }
    printf("d\n");

    p_im_state->E = compute_E_in(p_im_state->S, p_im_state->h_effective,
                                 p_in_model->h, p_in_model->num_nodes);
    printf("e\n");
}

void freeIMState(IMState* p_im_state, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params)
{

    free( p_im_state->S );
    free( p_im_state->h_effective );

    return;
}

void initSAWData(SAWData* p_saw_data, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params)
{
    // Also does the allocations
    int curr_node;
    double P;

    p_saw_data->p_down_heap = CreateSAWHeap( p_in_model->num_nodes );
    p_saw_data->p_up_heap = CreateSAWHeap( p_in_model->num_nodes );
    p_saw_data->sigma = (int*)malloc( p_in_model->num_nodes*sizeof(int) );
    p_saw_data->rho = (int*)malloc( p_in_model->num_nodes*sizeof(int) );

    p_saw_data->n_initial = 0;

    for (curr_node=0; curr_node<p_in_model->num_nodes; curr_node++)
    { 
        if (p_sim_params->S_ref[curr_node] != p_sim_params->S_in[curr_node])
        { 

            P = getPFlipGivenChosen(p_sim_params->h_effective_in,
                                    p_sim_params->S_in,
                                    p_saw_data->gamma_down, curr_node);

            AddToSAWHeap( p_saw_data->p_down_heap, P, curr_node );

            (p_saw_data->n_initial)++;

        }
        else
        { 

            P = getPFlipGivenChosen(p_sim_params->h_effective_in,
                                    p_sim_params->S_in, p_saw_data->gamma_up,
                                    curr_node);

            AddToSAWHeap( p_saw_data->p_up_heap, P, curr_node );

        }
    }

    return;
}

void freeSAWData(SAWData* p_saw_data, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params)
{

    p_saw_data->p_down_heap = DestroySAWHeap( p_saw_data->p_down_heap );
    p_saw_data->p_up_heap = DestroySAWHeap( p_saw_data->p_up_heap );

    free( p_saw_data->sigma );
    free( p_saw_data->rho );

    return;
}

void proposeIMMove(SAWData* p_saw_data, IMState* p_im_state,
                   const IMInputModel* p_in_model,
                   const IMSimParams* p_sim_params)
{
    // Main function to do the heapy flip proposal.
    int i, I, flipped_node;
    int* sigma = p_saw_data->sigma;
    int* rho = p_saw_data->rho;
    SAWHeap* p_down_heap = p_saw_data->p_down_heap;
    SAWHeap* p_up_heap = p_saw_data->p_up_heap;
    double P;
    int num_nodes = p_in_model->num_nodes;
    int num_edges = p_in_model->num_edges;
    const GraphNodes* graph_nodes = p_in_model->graph_nodes;
    const int** graph_edges = p_in_model->graph_edges;
    const double* J = p_in_model->J;

    p_saw_data->n_current = p_saw_data->n_initial;
    p_saw_data->log_f_fwd = 0.0;
    p_saw_data->log_f_rev = 0.0;

    // Down SAW:
    for (i=0; i<p_saw_data->SAW_length; i++)
    {
        /* Sample an element from the candidates heap according to the
           weights. This is not a graph node but an index in the heap */
        I = sampleFromSAWHeap( p_down_heap );

        /* Get the corresponding graph node: */
        flipped_node = ( (p_down_heap->SAWHeapArray)[I] ).graph_node_index;

        p_saw_data->log_f_fwd += log((p_down_heap->SAWHeapArray[I]).weight);
        p_saw_data->log_f_fwd -=
            log((p_down_heap->SAWHeapArray[0]).sum_curr_and_below);

        sigma[i] = flipped_node;

        removeFromSAWHeap(p_down_heap, I);

        /* First, update the effective fields due to the flip (the tabu
           status of variables has no bearing on this ) */
        singleFlipHTranslate((const int*)p_im_state->S,
                             graph_nodes, graph_edges,
                             p_im_state->h_effective, J, num_nodes,
                             flipped_node);

        /* Then call the function that adjusts the influence of the flip
           on the heap elements (namely the nbrs of flipped_node) */
        updateHeapAfterFlip(p_down_heap, (const int*)p_im_state->S,
                            (const double*)p_im_state->h_effective,
                            graph_nodes, graph_edges,
                            flipped_node, p_saw_data->gamma_down);
        updateHeapAfterFlip(p_up_heap, (const int*)p_im_state->S,
                            (const double*)p_im_state->h_effective,
                            graph_nodes, graph_edges,
                            flipped_node, p_saw_data->gamma_up);

        /* Update the state and the energy: */
        p_im_state->E += 2 * p_im_state->h_effective[flipped_node] *
                         p_im_state->S[flipped_node];
        p_im_state->S[flipped_node] *= -1;

        // Add to the "up" heap:
        P = getPFlipGivenChosen((const double*)p_im_state->h_effective,
                                (const int*)p_im_state->S,
                                p_saw_data->gamma_up, flipped_node);
        AddToSAWHeap(p_up_heap, P, flipped_node);

        p_saw_data->log_f_rev += log(P);
        p_saw_data->log_f_rev -= 
            log((p_up_heap->SAWHeapArray[0]).sum_curr_and_below);

        // The "set" we're on:
        (p_saw_data->n_current)--;
    }

    // Up SAW:
    for (i=0; i<p_saw_data->SAW_length; i++)
    { 
        /* Sample an element from the candidates heap according to the
           weights. This is not a graph node but an index in the heap */
        I = sampleFromSAWHeap( p_up_heap );

        /* Get the corresponding graph node: */
        flipped_node = ( (p_up_heap->SAWHeapArray)[I] ).graph_node_index;

        p_saw_data->log_f_fwd += log((p_up_heap->SAWHeapArray[I]).weight);
        p_saw_data->log_f_fwd -=
            log((p_up_heap->SAWHeapArray[0]).sum_curr_and_below);

        rho[i] = flipped_node;

        removeFromSAWHeap( p_up_heap, I );

        /* First, update the effective fields due to the flip (the tabu
           status of variables has no bearing on this ) */
        singleFlipHTranslate((const int*)p_im_state->S,
                             graph_nodes, graph_edges, p_im_state->h_effective,
                             J, num_nodes, flipped_node);

        /* Then call the function that adjusts the influence of the flip
           on the heap elements (namely the nbrs of flipped_node) */
        updateHeapAfterFlip(p_up_heap, (const int*)p_im_state->S,
                            (const double*)p_im_state->h_effective,
                            graph_nodes, graph_edges, flipped_node,
                            p_saw_data->gamma_up);
        updateHeapAfterFlip(p_down_heap, (const int*)p_im_state->S,
                            (const double*)p_im_state->h_effective,
                            graph_nodes, graph_edges, flipped_node,
                            p_saw_data->gamma_down);

        /* Update the state and the energy: */
        p_im_state->E += 2 * p_im_state->h_effective[flipped_node] *
                         p_im_state->S[flipped_node];
        p_im_state->S[flipped_node] *= -1;

        // Add to the "down" heap:
        P = getPFlipGivenChosen((const double*)p_im_state->h_effective,
                                (const int*)p_im_state->S,
                                p_saw_data->gamma_down, flipped_node);
        AddToSAWHeap(p_down_heap, P, flipped_node);

        p_saw_data->log_f_rev += log(P);
        p_saw_data->log_f_rev -=
            log((p_down_heap->SAWHeapArray[0]).sum_curr_and_below);

        // The "set" we're on:
        (p_saw_data->n_current)++;
    }

    return;
}

void undoIMMove(SAWData* p_saw_data, IMState* p_im_state,
                const IMInputModel* p_in_model,
                const IMSimParams* p_sim_params)
{
    // In the event of a rejection, we need to restore the state/heaps.
    // This is simply the "mirror-image" of proposeIMMove with no
    // sampling, but using the actual sequence of the proposal.

    int i, I, flipped_node;
    int* sigma = p_saw_data->sigma;
    int* rho = p_saw_data->rho;
    SAWHeap* p_down_heap = p_saw_data->p_down_heap;
    SAWHeap* p_up_heap = p_saw_data->p_up_heap;
    double P;
    int num_nodes = p_in_model->num_nodes;
    int num_edges = p_in_model->num_edges;
    const GraphNodes* graph_nodes = p_in_model->graph_nodes;
    const int** graph_edges = p_in_model->graph_edges;
    const double* J = p_in_model->J;

    // Reverse the up SAW:
    for (i=p_saw_data->SAW_length-1; i>=0; i-- ) { 

        flipped_node = rho[i];

        // Remove it from the down-heap:
        I = (p_down_heap->pReverseNodeIndex)[flipped_node];
        removeFromSAWHeap( p_down_heap, I );

        /* First, update the effective fields due to the flip (the tabu
           status of variables has no bearing on this ) */
        singleFlipHTranslate( (const int*)p_im_state->S, graph_nodes, graph_edges, p_im_state->h_effective, J, num_nodes, flipped_node );

        /* Then call the function that adjusts the influence of the flip
           on the heap elements (namely the nbrs of flipped_node) */
        updateHeapAfterFlip
            (p_up_heap,
             (const int*)p_im_state->S, (const double*)p_im_state->h_effective,
             graph_nodes, graph_edges, flipped_node, p_saw_data->gamma_up);
        updateHeapAfterFlip
            (p_down_heap,
             (const int*)p_im_state->S, (const double*)p_im_state->h_effective,
             graph_nodes, graph_edges, flipped_node, p_saw_data->gamma_down);

        /* Update the state and the energy: */
        p_im_state->E += 2*p_im_state->h_effective[flipped_node]*p_im_state->S[flipped_node];
        p_im_state->S[flipped_node] *= -1;

        // Return flipped_node to up-heap:
        P = getPFlipGivenChosen
            ((const double*)p_im_state->h_effective, (const int*)p_im_state->S,
             p_saw_data->gamma_up, flipped_node);
        AddToSAWHeap( p_up_heap, P, flipped_node );

        // The "set" we're on:
        (p_saw_data->n_current)--;
    }

    // Reverse the down SAW:
    for (i=p_saw_data->SAW_length-1; i>=0; i-- )
    { 
        flipped_node = sigma[i];

        // Remove it from the up-heap:
        I = (p_up_heap->pReverseNodeIndex)[flipped_node];
        removeFromSAWHeap( p_up_heap, I );

        /* First, update the effective fields due to the flip (the tabu
           status of variables has no bearing on this ) */
        singleFlipHTranslate( (const int*)p_im_state->S, graph_nodes, graph_edges, p_im_state->h_effective, J, num_nodes, flipped_node );

        /* Then call the function that adjusts the influence of the flip
           on the heap elements (namely the nbrs of flipped_node) */
        updateHeapAfterFlip
            (p_up_heap,
             (const int*)p_im_state->S, (const double*)p_im_state->h_effective,
             graph_nodes, graph_edges, flipped_node, p_saw_data->gamma_up);
        updateHeapAfterFlip
            (p_down_heap,
             (const int*)p_im_state->S, (const double*)p_im_state->h_effective,
             graph_nodes, graph_edges, flipped_node, p_saw_data->gamma_down);

        /* Update the state and the energy: */
        p_im_state->E += 2*p_im_state->h_effective[flipped_node]*p_im_state->S[flipped_node];
        p_im_state->S[flipped_node] *= -1;

        // Return flipped_node to down-heap:
        P = getPFlipGivenChosen
            ((const double*)p_im_state->h_effective, (const int*)p_im_state->S,
             p_saw_data->gamma_down, flipped_node);
        AddToSAWHeap( p_down_heap, P, flipped_node );

        // The "set" we're on:
        (p_saw_data->n_current)++;
    }

    return;
}

static double getPFlipGivenChosen
    (const double* h_effective, const int* S, double beta, int I)
{
    /* 
       21.Nov.07: 

       OVERFLOW SAFE.

       Do NOT confuse this with the VECTOR function used in
       sampleSAWSequence.c 

       This function merely computes the quantity:

       exp( -beta*h_effective[I]*S[I] ) 
       -----------------------------------------
       2*cosh( beta*h_effective[I]*S[I] )

       in an overflow-safe manner.

       If a process is uniformly chosing among non-tabu variables, this
       is the prob that the chosen variable flips.

       This value is NOT defined for tabu variables! You should never be
       calling this function on them. (In such a case, the returned
       value would simply be as if the variable was not tabu.)

     */

    double P, sum, logD, E_max;

    /* Which is bigger, exp( -beta*h_eff[I]*S[I] ) or 
       exp( beta*h_eff[I]*S[I] ) */
    E_max = -1.0*beta*h_effective[I]*S[I];
    if ( E_max < 0 ) {
        E_max = -1*E_max;
    }

    /* Calculation of the log of the denominator: works by "pulling
       out" the largest value: */

    /* The first element: */
    sum = exp( -1.0*beta*h_effective[I]*S[I] - E_max );
    /* The next one: */
    sum += exp( beta*h_effective[I]*S[I] - E_max );
    /* Log of the denominator: */
    logD = E_max + log(sum);

    /* Finally: */
    P = exp( -1.0*beta*h_effective[I]*S[I] - logD );

    return( P );
}

static void updateHeapAfterFlip
    (pSAWHeap p_candidates_heap, const int* S, const double* h_effective,
     const GraphNodes* graph_nodes, const int** graph_edges, int curr_node,
     double beta)
{
    /* The tabu-aware function to update the heap due to the influence
       of variable curr_node flipping. 

       Does not matter whether the flip has actually been done or not
       as long as h_effective has been adjusted for the NEIGHBOURS. 

       This is because this function only operates on the neighbours of
       curr_node.

       For each neighbor of curr_node, if the nbr is not tabu (check if
       not in the heap using pReverseNodeIndex) then we modify its weight
       due to the new h_effective using the changeSAWHeapElementKey
       function.

       If nbr is tabu, do nothing with that nbr

     */

    /* Usual misc graph stuff */
    int i, curr_nbr, v1, v2, curr_edge, num_nbrs, heap_array_el;
    int* incident_edges;

    double P;

    incident_edges = ( graph_nodes[curr_node] ).incident_edge_list;
    num_nbrs = ( graph_nodes[curr_node] ).num_nbrs;

    for (i=0; i<num_nbrs; i++) {

        curr_edge = incident_edges[i];
        v1 = graph_edges[curr_edge][0];
        v2 = graph_edges[curr_edge][1];

        if ( v1==curr_node ) {
            /* v2 is the nbr of curr_node: */
            curr_nbr = v2;
        } else {
            /* v1 is the nbr: */
            curr_nbr = v1;
        }
        /* The element in the heap array corresponding to curr_nbr (-1 if
           not in) */
        heap_array_el = (p_candidates_heap->pReverseNodeIndex)[curr_nbr];

        /* If curr_nbr is not tabu: */
        if ( heap_array_el != -1 ) {
            /* Calculate the new weight of P, and change the key of the el
               corresponding to curr_nbr */

            P = getPFlipGivenChosen( h_effective, S, beta, curr_nbr );
            changeSAWHeapElementKey( p_candidates_heap, heap_array_el, P );
        }

    } /* end for i = 0:num_nbrs-1 */

    return;
}

static void singleFlipHTranslate
    (const int* S, const GraphNodes* graph_nodes, const int** graph_edges,
     double* h_effective, const double* J, int num_nodes, int curr_node)
{
    /* Calculate the effective fields due to flipping bit curr_node in
       state S. The flip has NOT actually been done! */

    /* Usual misc graph stuff */
    int i, curr_nbr, v1, v2, curr_edge, num_nbrs;
    int* incident_edges;

    incident_edges = ( graph_nodes[curr_node] ).incident_edge_list;
    num_nbrs = ( graph_nodes[curr_node] ).num_nbrs;

    for (i=0; i<num_nbrs; i++) {

        curr_edge = incident_edges[i];
        v1 = graph_edges[curr_edge][0];
        v2 = graph_edges[curr_edge][1];

        if ( v1==curr_node ) {
            /* v2 is the nbr of curr_node: */
            curr_nbr = v2;
        } else {
            /* v1 is the nbr: */
            curr_nbr = v1;
        }

        /* Note "-" because S[curr_node] has NOT been flipped! */
        h_effective[curr_nbr] = h_effective[curr_nbr] - 2*S[curr_node]*J[curr_edge];

    }

    return;
}

void free_input_model(IMInputModel *in_model)
{
    destroy_graphedges(in_model->graph_edges, in_model->num_edges);
    destroy_graphnodes(in_model->graph_nodes);

    in_model->graph_edges = NULL;
    in_model->graph_nodes = NULL;
}

