import tables as t
import utils
import matplotlib.pyplot as plt
import numpy as np

def get_sard_tests_stats(graph_name = 'ISG3D_UAI', policy_num = 1000, res_num = 0, max_lag = 500, hdf5_fn='adaptive_saw.h5'):
    node = "/%s/%s/p%d/r%d" % (graph_name, 'SARD_EGO', policy_num, res_num)

    
    h5file = t.openFile(hdf5_fn, mode = "a")
    group = h5file.getNode(node)

    E_samples = group.E_samples[:]
    rewards = group.rewards[:]
    
    nadaptations = group.nadaptations.read()
    nruns_per_param_update = group.nruns_per_param_update.read()

    acf = map(lambda i: utils.compute_acf(E_samples[nruns_per_param_update*nadaptations:], i), xrange(max_lag+1))

    h5file.close()
    
    return E_samples, acf, [rewards[nruns_per_param_update*n+1][0] for n in xrange(nadaptations)]

def get_sw_stats(matfile_name='../Result/ferro_frus_sw80000_nonzero_h_temp1.mat', E_sample_key = 'energy_vec', max_lag = 500, start = 20000):
    import scipy.io as sio
    matfile = sio.loadmat(matfile_name)
    acf = map(lambda i: utils.compute_acf(matfile['energy_vec'][start:], i), xrange(max_lag+1))
    
    return matfile[E_sample_key], acf
    
def get_gibbs_stats(np_file = '../Result/Regular_Gibbs_ferrof_80000.npy', max_lag = 500, start = 20000):
    E_samples = np.load(np_file)
    acf = map(lambda i: utils.compute_acf(E_samples[start:], i), xrange(max_lag+1))
    
    return E_samples, acf

def get_bgibbs_stats(np_file = '../Result/Regular_Gibbs_ferrof_80000.npy', max_lag = 500, start = 20000):
    energies = np.load(np_file)
    E_samples = energies[0]
    total_acf = np.zeros(max_lag+1)
    for energy in energies:
	total_acf += map(lambda i: utils.compute_acf(energy[start:], i), xrange(max_lag+1))
    
    average_acf = total_acf/float(len(energies))
    
    return E_samples, average_acf

def plot_tests(gibbs_specs, sw_specs, sard_specs, y_min=-5500, y_max=-5900, fig_num_base = 0, save = True, RBM_bgibbs_specs = None, subsample = False):
    fig_num_base += 1
    plt.figure(fig_num_base)
    
    E_samples_sw, acf_sw = get_sw_stats(**sw_specs)
    sw_acf_plt = plt.plot(acf_sw, '-.b', linewidth=2.0)

    E_samples_sard, acf_sard, rewards = get_sard_tests_stats(**sard_specs)
    sard_acf_plt = plt.plot(acf_sard, '-r', linewidth=2.0)
    
    E_samples_gibbs, acf_gibbs = get_gibbs_stats(**gibbs_specs)
    gibbs_acf_plt = plt.plot(acf_gibbs, '--g', linewidth=2.0)

    if RBM_bgibbs_specs is not None:
	E_samples_bgibbs, acf_bgibbs = get_bgibbs_stats(**RBM_bgibbs_specs)
	bgibbs_acf_plt = plt.plot(acf_bgibbs, 'k:', linewidth=2.0)

    plt.ylim(ymax=1.0)
    plt.ylim(ymin=0.0)
    plt.xlabel('Lag', fontsize=16)
    plt.ylabel('Auto-correlation', fontsize=16)
    
    if RBM_bgibbs_specs is not None:
	plt.figlegend( (sard_acf_plt, sw_acf_plt, gibbs_acf_plt, bgibbs_acf_plt),
           ('SARDONICS', 'Swendsen-Wang', 'Gibbs', 'Block Gibbs')
           ,'upper left' )
    else:
	plt.figlegend( (sard_acf_plt, sw_acf_plt, gibbs_acf_plt),
           ('SARDONICS', 'Swendsen-Wang', 'Gibbs')
           ,'upper left' )

    if save:
	plot_fn = '../plots/'+sard_specs['graph_name'] + '_acf.pdf'
	plt.savefig(plot_fn, format='pdf')
    print "Produced plot %s" % (plot_fn)
    
    figsize=(10,10); subplot_base_num = 300
    
    fig_num_base += 1
    plt.figure(fig_num_base, figsize=figsize)
    plt.subplot(subplot_base_num+11)
    
    if not subsample:
	E_samples_sard_used = E_samples_sard[len(E_samples_sard)-20000:]
    else:
	E_samples_sard_used = E_samples_sard[::5]
    
    sard_energy = plt.plot(E_samples_sard_used, 'b')
    plt.xlabel('SARDONICS', fontsize=16)
    plt.ylim(ymax=y_max)
    plt.ylim(ymin=y_min)
    
    
    if RBM_bgibbs_specs is not None:
	plt.subplot(subplot_base_num+12)
	
	if not subsample:
	    bgibbs_energy_used = E_samples_bgibbs[len(E_samples_gibbs)-20000:]
	else:
	    bgibbs_energy_used = E_samples_bgibbs[::5]
	
	bgibbs_energy = plt.plot(bgibbs_energy_used, 'b')
	plt.xlabel('Block Gibbs', fontsize=16)
	plt.ylim(ymax=y_max)
	plt.ylim(ymin=y_min)
	plt.xlim(xmax=20000)
    else:
	plt.subplot(subplot_base_num+12)
	
	if not subsample:
	    E_sample_sw_used = E_samples_sw[len(E_samples_sw)-20000:]
	else:
	    E_sample_sw_used = E_samples_sw[::5]
	
	sw_energy = plt.plot(E_sample_sw_used, 'b')
	plt.xlabel('Swendsen-Wang', fontsize=16)
	plt.ylim(ymax=y_max)
	plt.ylim(ymin=y_min)
    
    plt.subplot(subplot_base_num+13)
    
    if not subsample:
	E_samples_gibbs_used = E_samples_gibbs[len(E_samples_gibbs)-20000:]
    else:
	E_samples_gibbs_used = E_samples_gibbs[::5]
    
    gibbs_energy = plt.plot(E_samples_gibbs_used, 'b')
    plt.xlabel('Gibbs', fontsize=16)
    plt.ylim(ymax=y_max)
    plt.ylim(ymin=y_min)
    
    
    if save:
	plot_fn = '../plots/'+sard_specs['graph_name'] + '_trace.pdf'
	plt.savefig(plot_fn, format='pdf')
    print "Produced plot %s" % (plot_fn)
    
    if RBM_bgibbs_specs is not None:
	fig_num_base += 1
	plt.figure(fig_num_base)
	if subsample:
	    E_sample_sw_used = E_samples_sw[len(E_samples_sw)-20000:]
	else:
	    E_sample_sw_used = E_samples_sw[::5]
	
	sw_energy = plt.plot(E_sample_sw_used, 'b')
	plt.xlabel('Iteration', fontsize=16)
	plt.ylabel('Energy', fontsize=16)

	if save:
	    plot_fn = '../plots/'+sard_specs['graph_name'] + '_sw_trace.pdf'
	    plt.savefig(plot_fn, format='pdf')
    
    fig_num_base += 1
    plt.figure(fig_num_base)
    plt.plot(rewards, '.')
    plt.xlabel('Number of Adaptations', fontsize=16)
    plt.ylabel('Reward', fontsize=16)
    
    if save:
	plot_fn = '../plots/'+sard_specs['graph_name'] + '_rewards.pdf'
	plt.savefig(plot_fn, format='pdf')
    print "Produced plot %s" % (plot_fn)
    
    if not save:
	plt.show()
    
    return fig_num_base
    
if __name__ == '__main__':
    max_lag = 2
    fig_num_base = 0
    
    # 1.0 frustrated 2d ISG
    #ferrof_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferrof_80000_temp1.npy', 'max_lag' : max_lag,}
    #ferrof_sw_specs = {'matfile_name' : '../Result/ferro_frus_sw80000_nonzero_h_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferrof_sard_specs = {'graph_name' : 'Ferro2D_Frustrated', 'policy_num' : 1000, 'res_num' : 15, 'max_lag' : max_lag,}
    #y_min = -5900
    #y_max = -5500
    
    #fig_num_base = plot_tests(ferrof_gibbs_specs, ferrof_sw_specs, ferrof_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    # Ising 3D temperature 0.8
    #isg3d_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_isg3d_80000.npy', 'max_lag' : max_lag,}
    #isg3d_sw_specs = {'matfile_name' : '../Result/isg3d_sw80000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #isg3d_sard_specs = {'graph_name' : 'ISG3D_UAI', 'policy_num' : 1000, 'res_num' : 8, 'max_lag' : max_lag,}
    #y_min = -1350
    #y_max = -1050
    
    #fig_num_base = plot_tests(isg3d_gibbs_specs, isg3d_sw_specs, isg3d_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    # Ising 3D temperature 1
    #isg3d_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_isg3d_80000_temp1.npy', 'max_lag' : max_lag,}
    #isg3d_sw_specs = {'matfile_name' : '../Result/isg3d_UAI_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #isg3d_sard_specs = {'graph_name' : 'ISG3D_UAI', 'policy_num' : 1000, 'res_num' : 12, 'max_lag' : max_lag,}
    #y_min = -1350
    #y_max = -1150
    
    #fig_num_base = plot_tests(isg3d_gibbs_specs, isg3d_sw_specs, isg3d_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    # Ferro Temperature 2.27
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_80000.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/Ferro_UAI_sw80000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 9, 'max_lag' : max_lag,}
    #y_min = -5800
    #y_max = -3900
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    # Ferro Temperature 5
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_80000_temp5.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/ferro_UAI_SW80000_temp5.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 18, 'max_lag' : max_lag,}
    #y_min = -2200
    #y_max = -1000
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    # Ferro Temperature 5 100000 subsample trace
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_100000_temp5.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/ferro_sw100000_temp5.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 18, 'max_lag' : max_lag,}
    #y_min = -2200
    #y_max = -1000
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base, subsample = True)
    
    # Ferro Temperature 2.27 100000 subsample trace
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_100000_tempcrit.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/ferro_sw100000_crittemp.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 13, 'max_lag' : max_lag, 'hdf5_fn' :'../Result/adaptive_saw.h5'}
    #y_min = -5800
    #y_max = -3900
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base, subsample = True)
    
    # Ferro Temperature 1
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_80000_temp1.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/ferro_UAI_SW80000_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 15, 'max_lag' : max_lag,}
    #y_min = -7500
    #y_max = -5500
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    #Ferro Temperature 1 100000 subsample trace
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_100000_temp1.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/ferro_sw100000_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 15, 'max_lag' : max_lag,}
    #y_min = -7500
    #y_max = -5500
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base, subsample = True)
    
    
    # 0.8 frustrated 2d ISG
    #ferrof_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferrof_80000.npy', 'max_lag' : max_lag,}
    #ferrof_sw_specs = {'matfile_name' : '../Result/ferro_frus_temp08_nonzeroh_sw80000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferrof_sard_specs = {'graph_name' : 'Ferro2D_Frustrated', 'policy_num' : 1000, 'res_num' : 18, 'max_lag' : max_lag, }
    #y_min = -5700
    #y_max = -5200

    #fig_num_base = plot_tests(ferrof_gibbs_specs, ferrof_sw_specs, ferrof_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    #RBM_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_rbm_80000.npy', 'max_lag' : max_lag,}
    #RBM_sw_specs = {'matfile_name' : '../Result/rbm_sw80000_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #RBM_sard_specs = {'graph_name' : 'RBM_UAI', 'policy_num' : 1000, 'res_num' : 8, 'max_lag' : max_lag, 'hdf5_fn' :'../Result/adaptive_saw.h5'}
    #RBM_bgibbs_specs = {'np_file' : '../Result/Block_RBM_3.npy', 'max_lag' : max_lag,}
    #y_min = -5800
    #y_max = -5500
    
    #fig_num_base = plot_tests(RBM_gibbs_specs, RBM_sw_specs, RBM_sard_specs, y_min, y_max, fig_num_base = fig_num_base, RBM_bgibbs_specs = RBM_bgibbs_specs)

    
    # 0.8 Chimera
    ferrof_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_chi_100000.npy', 'max_lag' : max_lag,}
    ferrof_sw_specs = {'matfile_name' : '../Result/chimera_sw100000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    ferrof_sard_specs = {'graph_name' : 'Chimera_Graph', 'policy_num' : 1000, 'res_num' : 3, 'max_lag' : max_lag, }
    y_min = -150
    y_max = -50

    fig_num_base = plot_tests(ferrof_gibbs_specs, ferrof_sw_specs, ferrof_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    # =================================================================================================================
    # Old code
    # =================================================================================================================

    #ferrof_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferrof_80000_temp1.npy', 'max_lag' : max_lag,}
    #ferrof_sw_specs = {'matfile_name' : '../Result/ferro_frus_sw80000_nonzero_h_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferrof_sard_specs = {'graph_name' : 'Ferro2D_Frustrated', 'policy_num' : 1000, 'res_num' : 12, 'max_lag' : max_lag,}
    #y_min = -5900
    #y_max = -5500
    
    #fig_num_base = plot_tests(ferrof_gibbs_specs, ferrof_sw_specs, ferrof_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    
    #ferrof_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferrof_80000.npy', 'max_lag' : max_lag,}
    #ferrof_sw_specs = {'matfile_name' : '../Result/ferro_frus_temp08_nonzeroh_sw80000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferrof_sard_specs = {'graph_name' : 'Ferro2D_Frustrated', 'policy_num' : 1000, 'res_num' : 11, 'max_lag' : max_lag,}
    #y_min = -5700
    #y_max = -5200

    #fig_num_base = plot_tests(ferrof_gibbs_specs, ferrof_sw_specs, ferrof_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    #isg3d_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_isg3d_80000.npy', 'max_lag' : max_lag,}
    #isg3d_sw_specs = {'matfile_name' : '../Result/isg3d_sw80000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #isg3d_sard_specs = {'graph_name' : 'ISG3D_UAI', 'policy_num' : 1000, 'res_num' : 2, 'max_lag' : max_lag,}
    #y_min = -1350
    #y_max = -1050
    
    #fig_num_base = plot_tests(isg3d_gibbs_specs, isg3d_sw_specs, isg3d_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    
    #ferro_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_ferro_80000.npy', 'max_lag' : max_lag,}
    #ferro_sw_specs = {'matfile_name' : '../Result/Ferro_UAI_sw80000.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #ferro_sard_specs = {'graph_name' : 'Ferro2D_UAI', 'policy_num' : 1000, 'res_num' : 1, 'max_lag' : max_lag,}
    #y_min = -5800
    #y_max = -3900
    
    #fig_num_base = plot_tests(ferro_gibbs_specs, ferro_sw_specs, ferro_sard_specs, y_min, y_max, fig_num_base = fig_num_base)
    
    
    #RBM_gibbs_specs = {'np_file' : '../Result/Regular_Gibbs_rbm_80000.npy', 'max_lag' : max_lag,}
    #RBM_sw_specs = {'matfile_name' : '../Result/rbm_sw80000_temp1.mat', 'E_sample_key' : 'energy_vec', 'max_lag' : max_lag,}
    #RBM_sard_specs = {'graph_name' : 'RBM_UAI', 'policy_num' : 1000, 'res_num' : 8, 'max_lag' : max_lag, 'hdf5_fn' :'../Result/adaptive_saw.h5'}
    #RBM_bgibbs_specs = {'np_file' : '../Result/Block_RBM_3.npy', 'max_lag' : max_lag,}
    #y_min = -5800
    #y_max = -5500
    
    #fig_num_base = plot_tests(RBM_gibbs_specs, RBM_sw_specs, RBM_sard_specs, y_min, y_max, fig_num_base = fig_num_base, RBM_bgibbs_specs = RBM_bgibbs_specs)

