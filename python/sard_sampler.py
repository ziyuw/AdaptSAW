from __future__ import division
from numpy import *
import ctypes as c

import graphs_ctypes as grph_c

libsardsampler = None

def load_sard_sampler_lib():
    import ASAW_config
    
    global libsardsampler
    if libsardsampler == None:
        #libsardsampler = c.cdll.LoadLibrary("/home/zwang/project/AdaptSAW/python/libsardsampler_ctypes.so")
        libsardsampler = c.cdll.LoadLibrary(ASAW_config.get_lib_dir() + "libsardsampler_ctypes.so")
            

def make_output_buffers(nmoves, iters_per_move_max):
    E_samples = zeros([nmoves, 1], dtype=double)
    E_proposeds = zeros([nmoves, 1], dtype=double)
    log_f_fwds = zeros([nmoves, 1], dtype=double)
    log_f_revs = zeros([nmoves, 1], dtype=double)
    MH_ratios = zeros([nmoves, 1], dtype=double)
    move_types = zeros([nmoves, 1], dtype=int32)

    used_sigma_lengths = \
        zeros([nmoves*iters_per_move_max, 1], dtype=int32)
    used_rho_lengths = \
        zeros([nmoves*iters_per_move_max, 1], dtype=int32)
    used_iters_per_move = zeros([nmoves, 1], dtype=int32)

    used_gammas_his = zeros([nmoves, 1], dtype=double)
    used_gammas_los = zeros([nmoves, 1], dtype=double)

    used_P_LL = zeros([nmoves, 1], dtype=double)
    used_P_LH = zeros([nmoves, 1], dtype=double)
    used_P_HL = zeros([nmoves, 1], dtype=double)
    
    rewards = zeros([nmoves, 1], dtype=double)

    return E_samples, E_proposeds, log_f_fwds, log_f_revs, \
           MH_ratios, move_types, \
           used_sigma_lengths, used_rho_lengths, \
           used_iters_per_move, used_gammas_his, used_gammas_los, \
           used_P_LL, used_P_LH, used_P_HL, rewards

def make_ctype_for_buffer(buffer):
    if buffer.dtype == dtype('double'):
        return buffer.ctypes.data_as(c.POINTER(c.c_double))
    elif buffer.dtype == dtype('int32'):
        return buffer.ctypes.data_as(c.POINTER(c.c_int32))
    else:
        raise "Invalid buffer dtype"

def make_ctype_for_var(var):
    if array(var).dtype == dtype('int'):
        return c.c_int32(var)
    elif array(var).dtype == dtype('double'):
        return c.c_double(var)
    else:
        raise "Invalid buffer dtype"
    
def SARDRun_Unif(nbrs_list_i, incident_edge_list_i, edges_i,
                 J, h, S_in, h_effective_in,
                 beta_true, gamma_hi_bounds, gamma_lo_bounds,
                 iters_per_move_bounds,
                 max_SAW_length_bounds, min_SAW_length_bounds,
                 nmoves, nruns_per_param_update, return_samples, strategy, P_LL = 0.0, P_HL = 0.0):

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
    S_in_c, S_end_c = \
        map(grph_c.make_ctypes_for_node_activations, [S_in, S_end])

    beta_true_c, \
    gamma_hi_max_c, gamma_hi_min_c, gamma_lo_max_c, gamma_lo_min_c, \
    max_SAW_length_bound_max_c, max_SAW_length_bound_min_c, \
    min_SAW_length_bound_max_c, min_SAW_length_bound_min_c, \
    iters_per_move_max_c, iters_per_move_min_c, \
    nmoves_c, nruns_per_param_update_c, P_LL_c, P_HL_c \
    = map(make_ctype_for_var,
          [beta_true,
           gamma_hi_bounds[1], gamma_hi_bounds[0],
           gamma_lo_bounds[1], gamma_lo_bounds[0],
           max_SAW_length_bounds[1], max_SAW_length_bounds[0],
           min_SAW_length_bounds[1], min_SAW_length_bounds[0],
           iters_per_move_bounds[1], iters_per_move_bounds[0],
           nmoves, nruns_per_param_update, P_LL, P_HL])
    
    iters_per_move_max = iters_per_move_bounds[1]

    E_samples, E_proposeds, log_f_fwds, log_f_revs, \
    MH_ratios, move_types, used_sigma_lengths, used_rho_lengths, \
    used_iters_per_move, used_gammas_his, used_gammas_los, \
    used_P_LL, used_P_LH, used_P_HL, _ \
    = make_output_buffers(nmoves, iters_per_move_max)

    E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c, \
    MH_ratios_c, move_types_c, \
    used_sigma_lengths_c, used_rho_lengths_c, used_iters_per_move_c, \
    used_gammas_his_c, used_gammas_los_c, \
    used_P_LL_c, used_P_LH_c, used_P_HL_c \
    = map(make_ctype_for_buffer,
          [E_samples, E_proposeds, log_f_fwds, log_f_revs,
           MH_ratios, move_types,
           used_sigma_lengths, used_rho_lengths, used_iters_per_move,
           used_gammas_his, used_gammas_los,
           used_P_LL, used_P_LH, used_P_HL])

    samples = []

    return_sample_lmbd = lambda sample, dim:\
            ego_c.return_sample(samples, sample, dim)

    # Create the callbacks.
    return_sample_c = ego_c.RETURN_SAMPLE(return_sample_lmbd)

    libsardsampler.sampleSARD_unif_ext\
        (num_nodes_c, num_edges_c,
         nbrs_list_ptrs_c, nbrs_list_lens_c,
         incident_edge_list_ptrs_c, incident_edge_list_lens_c,
         edges_i_c,
         J_c, h_c, S_in_c, S_end_c, h_effective_in_c,
         beta_true_c,
         gamma_hi_max_c, gamma_hi_min_c,
         gamma_lo_max_c, gamma_lo_min_c,
         max_SAW_length_bound_max_c, max_SAW_length_bound_min_c,
         min_SAW_length_bound_max_c, min_SAW_length_bound_min_c,
         iters_per_move_min_c, iters_per_move_max_c,
         nmoves_c, nruns_per_param_update_c, return_sample_c,
         P_LL_c, P_HL_c, strategy,
         E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c,
         MH_ratios_c, move_types_c,
         used_sigma_lengths_c, used_rho_lengths_c,
         used_iters_per_move_c,
         used_gammas_his_c, used_gammas_los_c,
         used_P_LL_c, used_P_LH_c, used_P_HL_c)
    
    if return_samples:
	return samples
    else:
	return S_end, E_samples, E_proposeds, log_f_fwds, log_f_revs, \
	    MH_ratios, move_types, \
	    used_sigma_lengths, used_rho_lengths, used_iters_per_move, \
	    used_gammas_his, used_gammas_los, \
	    used_P_LL, used_P_LH, used_P_HL

def SARDRun_Expert():
    pass

import ego_ctypes as ego_c
from ego.gaussianprocess import GaussianProcess
from ego.gaussianprocess.kernel import GaussianKernel_ard

def SARDRun_EGO(nbrs_list_i, incident_edge_list_i, edges_i,
                J, h, S_in, h_effective_in,
                beta_true, gamma_hi_bounds, gamma_lo_bounds,
                iters_per_move_bounds,
                max_SAW_length_bounds, min_SAW_length_bounds,
                nmoves, nadaptations, nruns_per_param_update,
                policy_size, strategy, exp_temp, MH_rate, return_samples,
                S_reset = None, h_effective_reset = None):
    assert(not(S_reset is None and h_effective_reset is not None))

    nparams = 8
    nmoves = array([nmoves], dtype='int')
    nadaptations = array([nadaptations], dtype='int')
    nruns_per_param_update = \
        array([nruns_per_param_update], dtype='int')
    policy_size = array([policy_size], dtype='int')

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
    S_in_c, S_reset_c, S_end_c = \
        map(grph_c.make_ctypes_for_node_activations,
            [S_in, S_reset, S_end])

    beta_true_c, \
    gamma_hi_max_c, gamma_hi_min_c, gamma_lo_max_c, gamma_lo_min_c, \
    max_SAW_length_bound_max_c, max_SAW_length_bound_min_c, \
    min_SAW_length_bound_max_c, min_SAW_length_bound_min_c, \
    iters_per_move_max_c, iters_per_move_min_c, \
    nmoves_c, nadaptations_c, nruns_per_param_update_c \
    = map(make_ctype_for_var,
          [beta_true,
           gamma_hi_bounds[1], gamma_hi_bounds[0],
           gamma_lo_bounds[1], gamma_lo_bounds[0],
           max_SAW_length_bounds[1], max_SAW_length_bounds[0],
           min_SAW_length_bounds[1], min_SAW_length_bounds[0],
           iters_per_move_bounds[1], iters_per_move_bounds[0],
           nmoves, nadaptations, nruns_per_param_update])

    iters_per_move_max = iters_per_move_bounds[1]

    E_samples, E_proposeds, log_f_fwds, log_f_revs, \
    MH_ratios, move_types, used_sigma_lengths, used_rho_lengths, \
    used_iters_per_move, used_gammas_his, used_gammas_los, \
    used_P_LL, used_P_LH, used_P_HL, rewards \
    = make_output_buffers(nmoves, iters_per_move_max)

    E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c, \
    MH_ratios_c, move_types_c, \
    used_sigma_lengths_c, used_rho_lengths_c, used_iters_per_move_c, \
    used_gammas_his_c, used_gammas_los_c, \
    used_P_LL_c, used_P_LH_c, used_P_HL_c, rewards_c \
    = map(make_ctype_for_buffer,
          [E_samples, E_proposeds, log_f_fwds, log_f_revs,
           MH_ratios, move_types,
           used_sigma_lengths, used_rho_lengths, used_iters_per_move,
           used_gammas_his, used_gammas_los,
           used_P_LL, used_P_LH, used_P_HL, rewards])

    return_samples_c = make_ctype_for_var(return_samples)
    policy_size_c = make_ctype_for_var(policy_size)
    policy = zeros([policy_size * nparams, 1], dtype=double)
    policy_c = make_ctype_for_buffer(policy)

    kernel_hyperparms = \
        array(
            [iters_per_move_bounds[1] - iters_per_move_bounds[0] + 1,
             min_SAW_length_bounds[1] - min_SAW_length_bounds[0] + 1,
             max_SAW_length_bounds[1] - max_SAW_length_bounds[0] + 1,
             1.0, 1.0, 1.0,
             gamma_hi_bounds[1] - gamma_hi_bounds[0],
             gamma_lo_bounds[1] - gamma_lo_bounds[0]])

    kernel = GaussianKernel_ard(kernel_hyperparms)
    gp = GaussianProcess(kernel, noise=0.1)

    # NOTE: This method is no longer useful. Take it out!
    def check_constraint(nx_pt, bounds, coords_size):
	if nx_pt[1] <=  nx_pt[2]:
	    return False
	else:
	    return False

    # Samples to be returned
    samples = []

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
            num_samples, check_constraint, exp_temp, MH_rate)

    return_sample_lmbd = \
	lambda sample_x_ptr, coords_size:\
            ego_c.return_sample(samples, sample_x_ptr, coords_size)

    # Create the callbacks.
    maximize_EI_c = ego_c.ACQ_MAXIMIZER_FUNC(maximize_EI_lmbd)
    add_data_to_gp_c = ego_c.GP_ADD_DATA_FUNC(add_data_to_gp_lmbd)
    maximize_gp_c = ego_c.GP_MAXIMIZER_FUNC(maximize_gp_lmbd)
    gp_posterior_mean_c = ego_c.GP_POSTERIOR_MEAN_FUNC(gp_posterior_mean_lmbd)
    maximize_EI_local_c = ego_c.LOCAL_ACQ_MAXIMIZER_FUNC(maximize_EI_local_lmbd)
    sample_gp_posterior_mean_c = ego_c.SAMPLE_GP_POSTERIOR_MEAN_FUNC(sample_gp_posterior_mean_lmbd)
    return_sample_c = ego_c.RETURN_SAMPLE(return_sample_lmbd)
    
    # Make the call.
    libsardsampler.sampleSARD_ego_ext\
        (num_nodes_c, num_edges_c,
         nbrs_list_ptrs_c, nbrs_list_lens_c,
         incident_edge_list_ptrs_c, incident_edge_list_lens_c,
         edges_i_c,
         J_c, h_c, S_in_c, S_reset_c, S_end_c,
         h_effective_in_c, h_effective_reset_c,
         beta_true_c,
         gamma_hi_max_c, gamma_hi_min_c,
         gamma_lo_max_c, gamma_lo_min_c,
         max_SAW_length_bound_max_c, max_SAW_length_bound_min_c,
         min_SAW_length_bound_max_c, min_SAW_length_bound_min_c,
         iters_per_move_min_c, iters_per_move_max_c,
         nmoves_c, nadaptations_c, nruns_per_param_update_c,
         maximize_EI_c, add_data_to_gp_c,
         maximize_gp_c, gp_posterior_mean_c, maximize_EI_local_c,
         sample_gp_posterior_mean_c, return_sample_c, return_samples_c,
         strategy, E_samples_c, E_proposeds_c, log_f_fwds_c, log_f_revs_c,
         MH_ratios_c, move_types_c,
         used_sigma_lengths_c, used_rho_lengths_c,
         used_iters_per_move_c,
         used_gammas_his_c, used_gammas_los_c,
         used_P_LL_c, used_P_LH_c, used_P_HL_c,
         rewards_c, policy_c, policy_size_c)

    print len(samples)

    return S_end, E_samples, E_proposeds, log_f_fwds, log_f_revs, \
           MH_ratios, move_types, \
           used_sigma_lengths, used_rho_lengths, used_iters_per_move, \
           used_gammas_his, used_gammas_los, \
           used_P_LL, used_P_LH, used_P_HL, \
           rewards, gp, policy, samples

