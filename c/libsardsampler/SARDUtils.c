#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<stdbool.h>
#include<mt19937ar.h>
#include "GraphNodes.h"
#include "SARDOutData.h"
#include "SAWHeaps.h"
#include "SARD.h"
#include "SARDUtils.h"

#include "utils.h"

static double getPFlipGivenChosen(const double* h_effective,
                                  const int* S, double beta, int I);
static void updateHeapAfterFlip
    (pSAWHeap p_candidates_heap, const int* S, const double* h_effective,
     const GraphNodes* graph_nodes, const int** graph_edges,
     int curr_node, double beta);
static void singleFlipHTranslate
    (const int* S,
     const GraphNodes* graph_nodes, const int** graph_edges,
     double* h_effective, const double* J,
     int num_nodes, int curr_node);

void initSARDState(SARDState* p_sard_state, const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params)
{
    // Also does the allocations
    p_sard_state->S = (int*)malloc(p_in_model->num_nodes*sizeof(int));
    p_sard_state->h_effective = (double*)malloc(p_in_model->num_nodes*sizeof(double));

    for (int i=0; i<p_in_model->num_nodes; i++)
    { 
        p_sard_state->S[i] = p_sim_params->S_in[i];
        p_sard_state->h_effective[i] = p_sim_params->h_effective_in[i];
    }

    p_sard_state->E = p_sim_params->E_in;
}

void setSARDState(SARDState* p_sard_state,
                const SARDInputModel* p_in_model,
                const int *S_new, const double* h_effective_new)
{
    for (int i = 0; i < p_in_model->num_nodes; i++)
    { 
        p_sard_state->S[i] = S_new[i];
        p_sard_state->h_effective[i] = h_effective_new[i];
    }

    p_sard_state->E = compute_E_in(p_sard_state->S, p_sard_state->h_effective,
                                   p_in_model->h, p_in_model->num_nodes);
}

void freeSARDState(SARDState* p_sard_state, const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params)
{
    free(p_sard_state->S);
    free(p_sard_state->h_effective);
}

void initSAWData(SAWData* p_saw_data, const SARDInputModel* p_in_model,
                 const SARDSimParams* p_sim_params)
{
    p_saw_data->p_saw_heap_L = NULL;
    p_saw_data->p_saw_heap_H = NULL;

    // Also does the allocations

    p_saw_data->sigma_lengths = (int*)malloc(p_sim_params->iters_per_move_max*sizeof(int));
    p_saw_data->rho_lengths = (int*)malloc(p_sim_params->iters_per_move_max*sizeof(int));

    p_saw_data->sigmas = (int**)malloc(p_sim_params->iters_per_move_max*sizeof(int*));
    p_saw_data->rhos = (int**)malloc(p_sim_params->iters_per_move_max*sizeof(int*));

    for (int i=0; i<p_sim_params->iters_per_move_max; i++)
    {
        p_saw_data->sigmas[i] =
            (int*)malloc( (p_sim_params->max_SAW_length_bound_max + p_sim_params->min_SAW_length_bound_max)*sizeof(int));
        p_saw_data->rhos[i] =
            (int*)malloc((p_sim_params->max_SAW_length_bound_max + p_sim_params->min_SAW_length_bound_max)*sizeof(int));
    }
}

void freeSAWData(SAWData* p_saw_data, const SARDInputModel* p_in_model,
                 const SARDSimParams* p_sim_params)
{
//     printf("iter_per_move: %d\n", p_sim_params->iters_per_move);
    if (p_saw_data->p_saw_heap_L != NULL)
    {
	DestroySAWHeap(p_saw_data->p_saw_heap_L);
    }
    
//     printf("iter_per_move: %d\n", p_sim_params->iters_per_move);
    
    if (p_saw_data->p_saw_heap_H != NULL)
    {
        DestroySAWHeap(p_saw_data->p_saw_heap_H);
    }
    
//     printf("iter_per_move: %d\n", p_sim_params->iters_per_move);
    for (int i=0; i<p_sim_params->iters_per_move_max; i++)
    {
	free(p_saw_data->sigmas[i]);
        free(p_saw_data->rhos[i]);
    }

//     printf("iter_per_move: %d\n", p_sim_params->iters_per_move);
    free(p_saw_data->sigmas);
    free(p_saw_data->rhos);
    free(p_saw_data->sigma_lengths);
    free(p_saw_data->rho_lengths);
}

void resetSAWHeapsCurrent(SARDState* p_sard_state, SAWData* p_saw_data, 
			  const SARDInputModel* p_in_model,
			  const SARDSimParams* p_sim_params)
{
    if (p_saw_data->p_saw_heap_L != NULL)
    {
        DestroySAWHeap(p_saw_data->p_saw_heap_L);
    }

    if (p_saw_data->p_saw_heap_H != NULL)
    {
        DestroySAWHeap(p_saw_data->p_saw_heap_H);
    }

    p_saw_data->p_saw_heap_L = CreateSAWHeap(p_in_model->num_nodes);
    p_saw_data->p_saw_heap_H = CreateSAWHeap(p_in_model->num_nodes);

    for (int curr_node=0; curr_node<p_in_model->num_nodes; curr_node++)
    { 
        // All nodes are in the heap
        double P_L = getPFlipGivenChosen(p_sard_state->h_effective,
                                         p_sard_state->S,
                                         p_sim_params->gamma_lo, curr_node);
        AddToSAWHeap(p_saw_data->p_saw_heap_L, P_L, curr_node);

        double P_H = getPFlipGivenChosen(p_sard_state->h_effective,
                                         p_sard_state->S,
                                         p_sim_params->gamma_hi, curr_node);
        // P_H = 1.0 - P_H;
        AddToSAWHeap(p_saw_data->p_saw_heap_H, P_H, curr_node);
    }
}

void resetSAWHeaps(SAWData* p_saw_data, const SARDInputModel* p_in_model,
                   const SARDSimParams* p_sim_params)
{
    if (p_saw_data->p_saw_heap_L != NULL)
    {
        DestroySAWHeap(p_saw_data->p_saw_heap_L);
    }

    if (p_saw_data->p_saw_heap_H != NULL)
    {
        DestroySAWHeap(p_saw_data->p_saw_heap_H);
    }

    p_saw_data->p_saw_heap_L = CreateSAWHeap(p_in_model->num_nodes);
    p_saw_data->p_saw_heap_H = CreateSAWHeap(p_in_model->num_nodes);

    for (int curr_node=0; curr_node<p_in_model->num_nodes; curr_node++)
    { 
        // All nodes are in the heap
        double P_L = getPFlipGivenChosen(p_sim_params->h_effective_in,
                                         p_sim_params->S_in,
                                         p_sim_params->gamma_lo, curr_node);
        AddToSAWHeap(p_saw_data->p_saw_heap_L, P_L, curr_node);

        double P_H = getPFlipGivenChosen(p_sim_params->h_effective_in,
                                         p_sim_params->S_in,
                                         p_sim_params->gamma_hi, curr_node);
        // P_H = 1.0 - P_H;
        AddToSAWHeap(p_saw_data->p_saw_heap_H, P_H, curr_node);
    }
}

void proposeSAWMove(SAWData* p_saw_data, SARDState* p_sard_state,
                    const SARDInputModel* p_in_model,
                    const SARDSimParams* p_sim_params)
{
    // Main function to do the heapy flip proposal.

    // The SAW lengths are selected in here and the full set of
    // iterations is gone through.

    // All proposal probabilities are calculated, but not the mixture.

    const int num_nodes = p_in_model->num_nodes;
    const int num_edges = p_in_model->num_edges;
    const GraphNodes* graph_nodes = p_in_model->graph_nodes;
    const int** graph_edges = p_in_model->graph_edges;
    const double* J = p_in_model->J;

    SAWHeap* p_saw_heap_L = p_saw_data->p_saw_heap_L;
    SAWHeap* p_saw_heap_H = p_saw_data->p_saw_heap_H;

    p_saw_data->log_f_fwd_LL = 0.0;
    p_saw_data->log_f_fwd_LH = 0.0;
    p_saw_data->log_f_fwd_HL = 0.0;

    p_saw_data->log_f_rev_LL = 0.0;
    p_saw_data->log_f_rev_LH = 0.0;
    p_saw_data->log_f_rev_HL = 0.0;

    // SIGMA sampling step:
    for (int curr_iter=0; curr_iter<p_sim_params->iters_per_move; curr_iter++)
    {
        p_saw_data->sigma_lengths[curr_iter] = p_sim_params->SAW_length_min + (int)( genrand_real2()*(p_sim_params->SAW_length_max-p_sim_params->SAW_length_min+1) );

        int sigma_length = p_saw_data->sigma_lengths[curr_iter];
        int* sigma = p_saw_data->sigmas[curr_iter];

        for (int i=0; i<sigma_length; i++)
        { 
            int I_L,I_H, flipped_node;

            if ((p_saw_data->move_type==LL) || (p_saw_data->move_type==LH))
            {
                // Sample from the L heap:
                I_L = sampleFromSAWHeap( p_saw_heap_L );
                /* Get the corresponding graph node: */
                flipped_node = ( (p_saw_heap_L->SAWHeapArray)[I_L] ).graph_node_index;
                // Index in the H heap:
                I_H = p_saw_heap_H->pReverseNodeIndex[flipped_node];
            }
            else
            { 
                // Sample from the H heap:
                I_H = sampleFromSAWHeap( p_saw_heap_H );
                /* Get the corresponding graph node: */
                flipped_node = ( (p_saw_heap_H->SAWHeapArray)[I_H] ).graph_node_index;
                // Index in the L heap:
                I_L = p_saw_heap_L->pReverseNodeIndex[flipped_node];
            }

            // Update forward flippies:
            p_saw_data->log_f_fwd_LL += log((p_saw_heap_L->SAWHeapArray[I_L] ).weight);
            p_saw_data->log_f_fwd_LL -= log((p_saw_heap_L->SAWHeapArray[0] ).sum_curr_and_below);

            p_saw_data->log_f_fwd_LH += log((p_saw_heap_L->SAWHeapArray[I_L] ).weight);
            p_saw_data->log_f_fwd_LH -= log((p_saw_heap_L->SAWHeapArray[0] ).sum_curr_and_below);

            p_saw_data->log_f_fwd_HL += log((p_saw_heap_H->SAWHeapArray[I_H] ).weight);
            p_saw_data->log_f_fwd_HL -= log((p_saw_heap_H->SAWHeapArray[0] ).sum_curr_and_below);

            sigma[i] = flipped_node;

            // Remove from both heaps:
            removeFromSAWHeap(p_saw_heap_L, I_L);
            removeFromSAWHeap(p_saw_heap_H, I_H);

            /* Update the effective fields due to the flip (the tabu
               status of variables has no bearing on this ) */
            singleFlipHTranslate((const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node);

            /* Then call the function that adjusts the influence of the flip
               on the heap elements (namely the nbrs of flipped_node) */
            updateHeapAfterFlip(p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip(p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            /* Update the state and the energy: */
            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;      
        }

        // Restore all nodes to heaps (order irrelevant)
        for (int i=0; i<sigma_length; i++)
        { 
            int flipped_node = sigma[i];
            // L-heap
            double P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node);
            AddToSAWHeap(p_saw_heap_L, P, flipped_node);
            // H-heap
            P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node);
            AddToSAWHeap(p_saw_heap_H, P, flipped_node);
        }

        // Calculation of log_f_rev. Dumb way by modifying the heap, though
        // in principle only h_effective is needed.
        for (int i=sigma_length-1; i>=0; i--)
        { 
            int flipped_node = sigma[i];
            int I_L = p_saw_heap_L->pReverseNodeIndex[flipped_node];
            int I_H = p_saw_heap_H->pReverseNodeIndex[flipped_node];

            // 5.Dec.10: consensus is that these are correct. For reference,
            // imagine that there was only one move type!
            p_saw_data->log_f_rev_LL += log((p_saw_heap_L->SAWHeapArray[I_L] ).weight);
            p_saw_data->log_f_rev_LL -= log((p_saw_heap_L->SAWHeapArray[0] ).sum_curr_and_below);

            p_saw_data->log_f_rev_LH += log((p_saw_heap_H->SAWHeapArray[I_H] ).weight);
            p_saw_data->log_f_rev_LH -= log((p_saw_heap_H->SAWHeapArray[0] ).sum_curr_and_below);

            p_saw_data->log_f_rev_HL += log((p_saw_heap_L->SAWHeapArray[I_L] ).weight);
            p_saw_data->log_f_rev_HL -= log((p_saw_heap_L->SAWHeapArray[0] ).sum_curr_and_below);

            removeFromSAWHeap(p_saw_heap_L, I_L);
            removeFromSAWHeap(p_saw_heap_H, I_H);

            singleFlipHTranslate((const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node);

            updateHeapAfterFlip(p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip(p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;
        }

        // Once again restore all nodes to heaps (order irrelevant)
        for (int i=0; i<sigma_length; i++)
        { 
            int flipped_node = sigma[i];
            // L-heap
            double P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node);
            AddToSAWHeap(p_saw_heap_L, P, flipped_node);
            // H-heap
            P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node);
            AddToSAWHeap(p_saw_heap_H, P, flipped_node);
        }

        // Bring the heaps back up. In this case no need to remove; can just
        // change the key of the flipped node.
        for (int i=0; i<sigma_length; i++)
        { 
            int flipped_node = sigma[i];
            singleFlipHTranslate((const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node);

            updateHeapAfterFlip(p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip(p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;

            int I_L= p_saw_heap_L->pReverseNodeIndex[flipped_node];
            double P_L = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node);
            changeSAWHeapElementKey(p_saw_heap_L, I_L, P_L);

            int I_H= p_saw_heap_H->pReverseNodeIndex[flipped_node];
            double P_H = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node);
            changeSAWHeapElementKey(p_saw_heap_H, I_H, P_H);
        }
    }

    // RHO sampling step:
    for (int curr_iter=0; curr_iter<p_sim_params->iters_per_move; curr_iter++)
    {
        p_saw_data->rho_lengths[curr_iter] = p_sim_params->SAW_length_min + (int)( genrand_real2()*(p_sim_params->SAW_length_max-p_sim_params->SAW_length_min+1) );

        int rho_length = p_saw_data->rho_lengths[curr_iter];
        int* rho = p_saw_data->rhos[curr_iter];

        for (int i=0; i<rho_length; i++)
        { 
            int I_L,I_H, flipped_node;

            if ((p_saw_data->move_type==LL) || (p_saw_data->move_type==HL))
            {
                // Sample from the L heap:
                I_L = sampleFromSAWHeap( p_saw_heap_L );
                /* Get the corresponding graph node: */
                flipped_node = ( (p_saw_heap_L->SAWHeapArray)[I_L] ).graph_node_index;
                // Index in the H heap:
                I_H = p_saw_heap_H->pReverseNodeIndex[flipped_node];
            }
            else
            { 
                // Sample from the H heap:
                I_H = sampleFromSAWHeap( p_saw_heap_H );
                /* Get the corresponding graph node: */
                flipped_node = ( (p_saw_heap_H->SAWHeapArray)[I_H] ).graph_node_index;
                // Index in the L heap:
                I_L = p_saw_heap_L->pReverseNodeIndex[flipped_node];
            }

            // Update forward flippies:
            p_saw_data->log_f_fwd_LL += log((p_saw_heap_L->SAWHeapArray[I_L] ).weight);
            p_saw_data->log_f_fwd_LL -= log((p_saw_heap_L->SAWHeapArray[0] ).sum_curr_and_below);

            p_saw_data->log_f_fwd_LH += log((p_saw_heap_H->SAWHeapArray[I_H] ).weight);
            p_saw_data->log_f_fwd_LH -= log((p_saw_heap_H->SAWHeapArray[0] ).sum_curr_and_below);

            p_saw_data->log_f_fwd_HL += log((p_saw_heap_L->SAWHeapArray[I_L] ).weight);
            p_saw_data->log_f_fwd_HL -= log((p_saw_heap_L->SAWHeapArray[0] ).sum_curr_and_below);

            rho[i] = flipped_node;

            // Remove from both heaps:
            removeFromSAWHeap(p_saw_heap_L, I_L);
            removeFromSAWHeap(p_saw_heap_H, I_H);

            /* Update the effective fields due to the flip (the tabu
               status of variables has no bearing on this ) */
            singleFlipHTranslate((const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node);

            /* Then call the function that adjusts the influence of the flip
               on the heap elements (namely the nbrs of flipped_node) */
            updateHeapAfterFlip(p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip(p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            /* Update the state and the energy: */
            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;      
        }
        // Restore all nodes to heaps (order irrelevant)
        for (int i=0; i<rho_length; i++)
        { 
            int flipped_node = rho[i];
            // L-heap
            double P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node);
            AddToSAWHeap(p_saw_heap_L, P, flipped_node);
            // H-heap
            P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node);
            AddToSAWHeap(p_saw_heap_H, P, flipped_node);
        }
        // Calculation of log_f_rev. Dumb way by modifying the heap, though
        // in principle only h_effective is needed.  
        for (int i=rho_length-1; i>=0; i--)
        { 
            int flipped_node = rho[i];
            int I_L = p_saw_heap_L->pReverseNodeIndex[flipped_node];
            int I_H = p_saw_heap_H->pReverseNodeIndex[flipped_node];

            // 5.Dec.10: consensus is that these are correct. For reference,
            // imagine that there was only one move type!
            p_saw_data->log_f_rev_LL += log((p_saw_heap_L->SAWHeapArray[I_L]).weight);
            p_saw_data->log_f_rev_LL -= log((p_saw_heap_L->SAWHeapArray[0]).sum_curr_and_below);

            p_saw_data->log_f_rev_LH += log((p_saw_heap_L->SAWHeapArray[I_L]).weight);
            p_saw_data->log_f_rev_LH -= log((p_saw_heap_L->SAWHeapArray[0]).sum_curr_and_below);

            p_saw_data->log_f_rev_HL += log((p_saw_heap_H->SAWHeapArray[I_H]).weight);
            p_saw_data->log_f_rev_HL -= log((p_saw_heap_H->SAWHeapArray[0]).sum_curr_and_below);

            removeFromSAWHeap(p_saw_heap_L, I_L);
            removeFromSAWHeap(p_saw_heap_H, I_H);

            singleFlipHTranslate((const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node);

            updateHeapAfterFlip(p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip(p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;
        }

        // Once again restore all nodes to heaps (order irrelevant)
        for (int i=0; i<rho_length; i++)
        { 
            int flipped_node = rho[i];
            // L-heap
            double P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node);
            AddToSAWHeap(p_saw_heap_L, P, flipped_node);
            // H-heap
            P = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node);
            AddToSAWHeap(p_saw_heap_H, P, flipped_node);
        }

        // Bring the heaps back up. In this case no need to remove; can just
        // change the key of the flipped node.
        for (int i=0; i<rho_length; i++)
        { 
            int flipped_node = rho[i];
            singleFlipHTranslate((const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node);

            updateHeapAfterFlip(p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip(p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;

            int I_L= p_saw_heap_L->pReverseNodeIndex[flipped_node];
            double P_L = getPFlipGivenChosen((const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node);
            changeSAWHeapElementKey(p_saw_heap_L, I_L, P_L);

            int I_H= p_saw_heap_H->pReverseNodeIndex[flipped_node];
            double P_H = getPFlipGivenChosen( (const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node );
            changeSAWHeapElementKey(p_saw_heap_H, I_H, P_H);
        }
    }
}


void undoSAWMove(SAWData* p_saw_data, SARDState* p_sard_state,
                 const SARDInputModel* p_in_model,
                 const SARDSimParams* p_sim_params)
{
    // Mirror image of propose. Restores state/heaps in the event of a rejection

    const int num_nodes = p_in_model->num_nodes;
    const int num_edges = p_in_model->num_edges;
    const GraphNodes* graph_nodes = p_in_model->graph_nodes;
    const int** graph_edges = p_in_model->graph_edges;
    const double* J = p_in_model->J;

    SAWHeap* p_saw_heap_L = p_saw_data->p_saw_heap_L;
    SAWHeap* p_saw_heap_H = p_saw_data->p_saw_heap_H;

    // Start with rho:
    for ( int curr_iter=p_sim_params->iters_per_move-1; curr_iter>=0; curr_iter-- )
    {
        int rho_length = p_saw_data->rho_lengths[curr_iter];
        const int* rho = p_saw_data->rhos[curr_iter];

        for ( int i=rho_length-1; i>=0; i-- ) { 

            int flipped_node = rho[i];

            /* Update the effective fields due to the flip (the tabu status of
               variables has no bearing on this ) */
            singleFlipHTranslate( (const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node );

            /* Then call the function that adjusts the influence of the flip
               on the heap elements (namely the nbrs of flipped_node) */
            updateHeapAfterFlip( p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo);
            updateHeapAfterFlip( p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi);

            /* Update the state and the energy: */
            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;

            // Change the keys (no need to remove/add)
            int I_L = (p_saw_heap_L->pReverseNodeIndex)[flipped_node];
            double P_L = getPFlipGivenChosen( (const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node );
            changeSAWHeapElementKey( p_saw_heap_L, I_L, P_L );

            int I_H = (p_saw_heap_H->pReverseNodeIndex)[flipped_node];
            double P_H = getPFlipGivenChosen( (const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node );
            changeSAWHeapElementKey( p_saw_heap_H, I_H, P_H );

        }
    }

    // Now sigma:
    for ( int curr_iter=p_sim_params->iters_per_move-1; curr_iter>=0; curr_iter-- )
    {
        int sigma_length = p_saw_data->sigma_lengths[curr_iter];
        const int* sigma = p_saw_data->sigmas[curr_iter];

        for ( int i=sigma_length-1; i>=0; i-- ) { 

            int flipped_node = sigma[i];

            /* Update the effective fields due to the flip (the tabu status of
               variables has no bearing on this ) */
            singleFlipHTranslate( (const int*)p_sard_state->S, graph_nodes, graph_edges, p_sard_state->h_effective, J, num_nodes, flipped_node );

            /* Then call the function that adjusts the influence of the flip
               on the heap elements (namely the nbrs of flipped_node) */
            updateHeapAfterFlip( p_saw_heap_L, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_lo );
            updateHeapAfterFlip( p_saw_heap_H, (const int*)p_sard_state->S, (const double*)p_sard_state->h_effective, graph_nodes, graph_edges, flipped_node, p_sim_params->gamma_hi );

            /* Update the state and the energy: */
            p_sard_state->E += 2*p_sard_state->h_effective[flipped_node]*p_sard_state->S[flipped_node];
            p_sard_state->S[flipped_node] *= -1;

            // Change the keys (no need to remove/add)
            int I_L = (p_saw_heap_L->pReverseNodeIndex)[flipped_node];
            double P_L = getPFlipGivenChosen( (const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_lo, flipped_node );
            changeSAWHeapElementKey( p_saw_heap_L, I_L, P_L );

            int I_H = (p_saw_heap_H->pReverseNodeIndex)[flipped_node];
            double P_H = getPFlipGivenChosen( (const double*)p_sard_state->h_effective, (const int*)p_sard_state->S, p_sim_params->gamma_hi, flipped_node );
            changeSAWHeapElementKey( p_saw_heap_H, I_H, P_H );

        }
    }
}

static double getPFlipGivenChosen(const double* h_effective, const int* S,
                                  double beta, int I)
{
    /* 
       21.Nov.07: 

       OVERFLOW SAFE.

       Do NOT confuse this with the VECTOR function used in
       sampleSAWSequence.c 

       In an overflow-safe manner:

       computes the quantity:

       exp( -beta*h_effective[I]*S[I] ) 
       -----------------------------------------
       2*cosh( beta*h_effective[I]*S[I] )

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

    P = exp( -1.0*beta*h_effective[I]*S[I] - logD );

    return( P );
}


static void updateHeapAfterFlip
    (pSAWHeap p_candidates_heap, const int* S,const double* h_effective,
     const GraphNodes* graph_nodes, const int** graph_edges,
     int curr_node, double beta)
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

  for (i=0; i<num_nbrs; i++)
  {
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
}

void evalLogFs(double* p_log_f_fwd, double* p_log_f_rev,
                      const SAWData* p_saw_data,
                      const SARDSimParams* p_sim_params)
{
    // Overflow-safe.

    double log_P_LL = log(p_sim_params->P_LL),
           log_P_LH = log(p_sim_params->P_LH),
           log_P_HL = log(p_sim_params->P_HL);
    double A,B,C;

    // First fwd:
    A = log_P_LL + p_saw_data->log_f_fwd_LL;
    B = log_P_LH + p_saw_data->log_f_fwd_LH;
    C = log_P_HL + p_saw_data->log_f_fwd_HL;

    double max = A;
    if ( B > max ) { 
        max = B;
    }
    if ( C > max ) { 
        max = C;
    }

    *p_log_f_fwd = max + log( exp(A-max) + exp(B-max) + exp(C-max) );

    // Now revs:
    A = log_P_LL + p_saw_data->log_f_rev_LL;
    B = log_P_LH + p_saw_data->log_f_rev_LH;
    C = log_P_HL + p_saw_data->log_f_rev_HL;

    max = A;
    if ( B > max ) { 
        max = B;
    }
    if ( C > max ) { 
        max = C;
    }

    *p_log_f_rev = max + log( exp(A-max) + exp(B-max) + exp(C-max) );
}

