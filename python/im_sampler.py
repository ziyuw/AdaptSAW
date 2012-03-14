from __future__ import division
from numpy import *
import ctypes as c

import graphs_ctypes as grph_c

libimsampler = None

def load_im_sampler_lib():
    global libimsampler
    if libimsampler == None:
	import ASAW_config
        libimsampler = \
	    c.cdll.LoadLibrary(ASAW_config.get_lib_dir() + "libimsampler_ctypes.so")

def make_output_buffers(nmoves, nactions):
    E_samples = zeros([nmoves, 1], dtype=double)
    E_proposeds = zeros([nmoves, 1], dtype=double)
    log_f_fwds = zeros([nmoves, 1], dtype=double)
    log_f_revs = zeros([nmoves, 1], dtype=double)
    MH_ratios = zeros([nmoves, 1], dtype=double)
    SAW_lengths = zeros([nmoves, 1], dtype=int32)
    used_gammas = zeros([nmoves, 1], dtype=double)
    rewards = zeros([nmoves, 1], dtype=double)

    return E_samples, E_proposeds, log_f_fwds, log_f_revs, \
           MH_ratios, SAW_lengths, used_gammas, rewards

def make_ctypes_for_output_buffers\
        (E_samples, E_proposeds, log_f_fwds, log_f_revs,
         MH_ratios, SAW_lengths, used_gammas, rewards):
    E_samples_c = E_samples.ctypes.data_as(c.POINTER(c.c_double))
    E_proposeds_c = E_proposeds.ctypes.data_as(c.POINTER(c.c_double))
    log_f_fwds_c = log_f_fwds.ctypes.data_as(c.POINTER(c.c_double))
    log_f_revs_c = log_f_revs.ctypes.data_as(c.POINTER(c.c_double))
    MH_ratios_c = MH_ratios.ctypes.data_as(c.POINTER(c.c_double))
    SAW_lengths_c = SAW_lengths.ctypes.data_as(c.POINTER(c.c_int32))
    used_gammas_c = used_gammas.ctypes.data_as(c.POINTER(c.c_double))
    rewards_c = rewards.ctypes.data_as(c.POINTER(c.c_double))

    return E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c, \
           MH_ratios_c, SAW_lengths_c, used_gammas_c, rewards_c

def make_ctypes_for_simulation_parameters\
        (beta_true, nmoves, nruns_per_param_update, gammas,
         SAW_length_min, SAW_length_max):
    beta_true_c = c.c_double(beta_true)

    gammas_c = gammas.ctypes.data_as(c.POINTER(c.c_double))
    len_gammas_c = c.c_int32(len(gammas))

    SAW_length_min_c = c.c_int32(SAW_length_min)
    SAW_length_max_c = c.c_int32(SAW_length_max)

    nmoves_c = c.c_int32(nmoves)
    nruns_per_param_update_c = c.c_int32(nruns_per_param_update)

    return beta_true_c, nmoves_c, nruns_per_param_update_c, \
           gammas_c, len_gammas_c, \
           SAW_length_min_c, SAW_length_max_c

def IMRun_Unif(nbrs_list_i, incident_edge_list_i, edges_i,
               J, h, S_ref, S_in, h_effective_in,
               beta_true, gammas, SAW_length_min, SAW_length_max,
               nmoves, nruns_per_param_update):
    num_nodes_c, num_edges_c,\
    edges_i_c, \
    nbrs_list_ptrs, nbrs_list_ptrs_c, \
    nbrs_list_lens, nbrs_list_lens_c, \
    incident_edge_list_ptrs, incident_edge_list_ptrs_c, \
    incident_edge_list_lens, incident_edge_list_lens_c \
    = grph_c.make_ctypes_for_graph_repr \
      (nbrs_list_i, incident_edge_list_i, edges_i)

    J_c, h_c, h_effective_in_c = \
        map(grph_c.make_ctypes_for_graph_weights,
            [J, h, h_effective_in])

    S_end = zeros(S_in.shape, dtype=int32)
    S_ref_c, S_in_c, S_end_c = \
        map(grph_c.make_ctypes_for_node_activations,
            [S_ref, S_in, S_end])

    beta_true_c, nmoves_c, nruns_per_param_update_c, \
    gammas_c, len_gammas_c, \
    SAW_length_min_c, SAW_length_max_c \
    = make_ctypes_for_simulation_parameters\
      (beta_true, nmoves, nruns_per_param_update, \
       gammas, SAW_length_min, SAW_length_max)

    E_samples, E_proposeds, log_f_fwds, log_f_revs, \
    MH_ratios, used_gammas = \
       [zeros([nmoves, 1], dtype=double) for i in xrange(6)]
    
    SAW_lengths = zeros([nmoves, 1], dtype=int32)

    E_samples_c = E_samples.ctypes.data_as(c.POINTER(c.c_double))
    E_proposeds_c = E_proposeds.ctypes.data_as(c.POINTER(c.c_double))
    log_f_fwds_c = log_f_fwds.ctypes.data_as(c.POINTER(c.c_double))
    log_f_revs_c = log_f_revs.ctypes.data_as(c.POINTER(c.c_double))
    MH_ratios_c = MH_ratios.ctypes.data_as(c.POINTER(c.c_double))
    used_gammas_c = used_gammas.ctypes.data_as(c.POINTER(c.c_double))

    SAW_lengths_c = SAW_lengths.ctypes.data_as(c.POINTER(c.c_int32))

    libimsampler.sampleIM_unif_ext\
        (num_nodes_c, num_edges_c,
         nbrs_list_ptrs_c, nbrs_list_lens_c,
         incident_edge_list_ptrs_c, incident_edge_list_lens_c,
         edges_i_c,
         J_c, h_c, S_ref_c, S_in_c, S_end_c, h_effective_in_c,
         beta_true_c,
         len_gammas_c, gammas_c,
         SAW_length_min_c, SAW_length_max_c,
         nmoves_c, nruns_per_param_update_c,
         E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c,
         MH_ratios_c, SAW_lengths_c, used_gammas_c)

    # Return the answers.
    return S_end, E_samples, E_proposeds, log_f_fwds, log_f_revs, \
           MH_ratios, SAW_lengths, used_gammas

def Kawasaki_Sampler(nbrs_list_i, incident_edge_list_i, edges_i,
                     J, h, S_ref, S_in, h_effective_in,
                     beta_true, nmoves):
    gammas = resize(array(float(0)), (1,))
    SAW_length = 1
    return IMRun_Unif(nbrs_list_i, incident_edge_list_i, edges_i,
                      J, h, S_ref, S_in, h_effective_in,
                      beta_true, gammas, SAW_length, SAW_length,
                      nmoves, nmoves)

def IMRun_Expert(nbrs_list_i, incident_edge_list_i, edges_i,
                 J, h, S_ref, S_in, h_effective_in,
                 beta_true, gamma, SAW_length_min, SAW_length_max,
                 nmoves, nruns_per_param_update):
    gammas = resize(array(float(gamma)), (1,))
    return IMRun_Unif(nbrs_list_i, incident_edge_list_i, edges_i,
                      J, h, S_ref, S_in, h_effective_in,
                      beta_true, gammas, SAW_length_min, SAW_length_max,
                      nmoves, nruns_per_param_update)

import ego_ctypes as ego_c
from ego.gaussianprocess import GaussianProcess
from ego.gaussianprocess.kernel import GaussianKernel_ard

def IMRun_EGO(nbrs_list_i, incident_edge_list_i, edges_i,
              J, h, S_ref, S_in, h_effective_in,
              beta_true, gammas, SAW_length_min, SAW_length_max,
              nmoves, nadaptations, nruns_per_param_update,
              policy_size, strategy, exp_temp, MH_rate, 
              S_reset = None, h_effective_reset = None):
    assert(not(S_reset is None and h_effective_reset is not None))

    print 'strategy', strategy

    if S_reset is None:
        S_reset = S_in
        h_effective_reset = h_effective_in
        
    num_nodes_c, num_edges_c,\
    edges_i_c, \
    nbrs_list_ptrs, nbrs_list_ptrs_c, \
    nbrs_list_lens, nbrs_list_lens_c, \
    incident_edge_list_ptrs, incident_edge_list_ptrs_c, \
    incident_edge_list_lens, incident_edge_list_lens_c \
    = grph_c.make_ctypes_for_graph_repr \
      (nbrs_list_i, incident_edge_list_i, edges_i)

    J_c, h_c, h_effective_in_c, h_effective_reset_c = \
        map(grph_c.make_ctypes_for_graph_weights,
            [J, h, h_effective_in, h_effective_reset])

    S_end = zeros(S_in.shape, dtype=int32)
    S_ref_c, S_in_c, S_reset_c, S_end_c = \
        map(grph_c.make_ctypes_for_node_activations,
            [S_ref, S_in, S_reset, S_end])

    beta_true_c, nmoves_c, nruns_per_param_update_c, \
    gammas_c, len_gammas_c, \
    SAW_length_min_c, SAW_length_max_c \
    = make_ctypes_for_simulation_parameters\
      (beta_true, nmoves, nruns_per_param_update, \
       gammas, SAW_length_min, SAW_length_max)

    nadaptations_c = c.c_int32(nadaptations)

    nactions = (SAW_length_max - SAW_length_min + 1) * len(gammas)

    E_samples, E_proposeds, log_f_fwds, log_f_revs, \
    MH_ratios, SAW_lengths, used_gammas, rewards \
    = make_output_buffers(nmoves, nactions)

    E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c, \
    MH_ratios_c, SAW_lengths_c, used_gammas_c, rewards_c \
    = make_ctypes_for_output_buffers\
        (E_samples, E_proposeds, log_f_fwds, log_f_revs,
         MH_ratios, SAW_lengths, used_gammas, rewards)
    
    policy_size_c = c.c_int32(policy_size)
    policy = zeros([policy_size*2, 1], dtype=double)
    policy_c = policy.ctypes.data_as(c.POINTER(c.c_double))
    
    kernel_hyperparms = \
        array([SAW_length_max - SAW_length_min + 1,
               gammas[-1] - gammas[0]]) * 0.1
    
    kernel = GaussianKernel_ard(kernel_hyperparms)
    gp = GaussianProcess(kernel, noise=0.1)

    # Bind the GP in closures. References to these closures need to be
    # kept, so that they only get garbage-collected when the C code
    # has finished running.
    maximize_EI_lmbd = \
        lambda bounds_ptr, max_x_ptr, x_size:\
            ego_c.maximizeEI_wrap(gp, bounds_ptr, max_x_ptr, x_size)

    add_data_to_gp_lmbd = \
        lambda x_ptr, x_size, y:\
            ego_c.add_data_to_gp(gp, x_ptr, x_size, y)

    maximize_gp_lmbd = \
        lambda x_ptr, x_size, y_ptr:\
            ego_c.maximize_gp(gp, x_ptr, x_size, y_ptr)

    gp_posterior_mean_lmbd = \
        lambda coords_ptr, coords_size:\
            ego_c.gp_posterior_mean(gp, coords_ptr, coords_size)

    maximize_EI_local_lmbd = \
        lambda bounds_ptr, max_x_ptr, x_size:\
            ego_c.maximizeEI_local_wrap(gp, bounds_ptr, max_x_ptr, x_size, strategy)
    
    sample_gp_posterior_mean_lmbd = \
	lambda coords_ptr, coords_size, bounds_ptr, sample_x_ptr, num_samples:\
            ego_c.sample_gp_posterior_mean(gp, coords_ptr, coords_size, bounds_ptr, sample_x_ptr, 
            num_samples, None, exp_temp, MH_rate)

    # Create the callbacks.
    maximize_EI_c = ego_c.ACQ_MAXIMIZER_FUNC(maximize_EI_lmbd)
    add_data_to_gp_c = ego_c.GP_ADD_DATA_FUNC(add_data_to_gp_lmbd)
    maximize_gp_c = ego_c.GP_MAXIMIZER_FUNC(maximize_gp_lmbd)
    gp_posterior_mean_c = ego_c.GP_POSTERIOR_MEAN_FUNC(gp_posterior_mean_lmbd)
    maximize_EI_local_c = ego_c.LOCAL_ACQ_MAXIMIZER_FUNC(maximize_EI_local_lmbd)
    sample_gp_posterior_mean_c = ego_c.SAMPLE_GP_POSTERIOR_MEAN_FUNC(sample_gp_posterior_mean_lmbd)

    # Make the call.
    libimsampler.sampleIM_ego_ext\
        (num_nodes_c, num_edges_c,
         nbrs_list_ptrs_c, nbrs_list_lens_c,
         incident_edge_list_ptrs_c, incident_edge_list_lens_c,
         edges_i_c,
         J_c, h_c,
         S_ref_c, S_in_c, S_reset_c, S_end_c,
         h_effective_in_c, h_effective_reset_c,
         beta_true_c,
         len_gammas_c, gammas_c,
         SAW_length_min_c, SAW_length_max_c,
         nmoves_c, nadaptations_c, nruns_per_param_update_c,
         maximize_EI_c, add_data_to_gp_c, maximize_gp_c, gp_posterior_mean_c,
         maximize_EI_local_c, sample_gp_posterior_mean_c,
         strategy, E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c,
         MH_ratios_c, SAW_lengths_c, used_gammas_c, rewards_c,
         policy_c, policy_size_c)

    gamma_min = gammas[0]
    gamma_max = gammas[-1]

    # Return the answers.
    return S_end, E_samples, E_proposeds, log_f_fwds, log_f_revs, \
           MH_ratios, SAW_lengths, used_gammas, rewards, gp, \
           gamma_min, gamma_max, policy

import ising
def IMRun_EGO_premixed(nbrs_list_i, incident_edge_list_i, edges_i,
                       J, h, S_ref, S_in, h_effective_in,
                       beta_true, gammas, SAW_length_min, SAW_length_max,
                       nmoves, nadaptations, nruns_per_param_update,
                       policy_size,
                       nadaptations_premixing = None,
                       nruns_per_param_update_premixing = None):
    strategy = 1
			   
    if nadaptations_premixing is None:
        nadaptations_premixing = nadaptations
    if nruns_per_param_update_premixing is None:
        nruns_per_param_update_premixing = nruns_per_param_update
    
    S_start = \
        IM_chainmixer_EGO(nbrs_list_i, incident_edge_list_i, edges_i,
                          J, h, S_ref, S_in, h_effective_in,
                          beta_true, gammas, SAW_length_min, SAW_length_max,
                          nadaptations_premixing,
                          nruns_per_param_update_premixing, strategy)

    h_effective_start = \
        ising.get_effective_fields(edges_i, incident_edge_list_i,
                                   J, h, S_start)

    return IMRun_EGO(nbrs_list_i, incident_edge_list_i, edges_i,
                     J, h, S_ref, S_start, h_effective_start,
                     beta_true, gammas, SAW_length_min, SAW_length_max,
                     nmoves, nadaptations, nruns_per_param_update,
                     policy_size, strategy,
                     S_in, h_effective_in)

# Use EGO and IM to find a low-energy state that is far from the
# initial state. This is done by running the IM_EGO sampler for a
# short adaptation phase and then using the state that the sampler
# ended up in.
def IM_chainmixer_EGO(nbrs_list_i, incident_edge_list_i, edges_i,
                      J, h, S_ref, S_in, h_effective_in,
                      beta_true, gammas, SAW_length_min, SAW_length_max,
                      nadaptations, nruns_per_param_update, strategy):
    throwaway_policy_size = 100
    nmoves = nadaptations * nruns_per_param_update
    S_start = IMRun_EGO(nbrs_list_i, incident_edge_list_i, edges_i,
                        J, h, S_ref, S_in, h_effective_in,
                        beta_true, gammas, SAW_length_min, SAW_length_max,
                        nmoves, nadaptations, nruns_per_param_update,
                        throwaway_policy_size, strategy)[0]
    return S_start
