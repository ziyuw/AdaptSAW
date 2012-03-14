from __future__ import division

import utils
utils.load_matplotlib_backend()
        
import graphs as g
import ising
import im_sampler as im
import sard_sampler as sard

from numpy import *

def get_num_adapted_moves(nruns_per_param_update, nadaptations):
    return nruns_per_param_update * nadaptations

import tables as t
hdf5_fn = 'adaptive_saw.h5'

def save_to_hdf5_file(loc, name, data):
    h5file = t.openFile(hdf5_fn, mode='a')

    try:
        node = h5file.getNode(loc, name)
        node._f_remove()
    except t.NoSuchNodeError as e:
        pass
    
    h5file.createArray(loc, name, data, '', None, True)
    h5file.flush()

    h5file.close()
    h5file = None

def create_next_run_group(h5file, loc):
    try:
        nruns = get_num_run_groups(h5file, loc)
    except t.NoSuchNodeError as e:
        nruns = 0

    run_name = 'r' + str(nruns)
    run_grp = h5file.createGroup(loc, run_name, '', None, True)

    return run_grp, nruns

def get_num_run_groups(h5file, loc):
    grp = h5file.getNode(loc)
    nruns = len(grp._v_groups)
    return nruns

class Sampler_Test(object):
    name = None
    hdf5_base_loc = None

    S_in = None
    nmoves = None
    
    res_num = None

    saved_res_fields = \
        set(['S_end','S_in',
             'E_samples', 'E_proposeds', 'log_f_fwds', 'log_f_revs',
             'MH_ratios', 'nmoves'])

    graph = None

    def run(self):
        pass

    def get_last_res_num(self, hdf5_fn='adaptive_saw.h5'):
	h5file = t.openFile(hdf5_fn, mode='a')
	res_num = get_num_run_groups(h5file, self.hdf5_base_loc)
	
	return res_num

    def get_hdf5_location(self, res_num):
        return self.hdf5_base_loc + "/r" + str(int(res_num))

    def load_results(self, res_num=None, hdf5_fn='adaptive_saw.h5',
                     force=False):
        h5file = t.openFile(hdf5_fn, mode='a')

        if res_num is None:
            res_num = get_num_run_groups(h5file, self.hdf5_base_loc) - 1
            
        run_grp_loc = self.get_hdf5_location(res_num)
        run_grp = h5file.getNode(run_grp_loc)
        
        map(lambda x: setattr(self, x, run_grp._v_leaves[x].read()),
            self.saved_res_fields)
        self.res_num = res_num

        h5file.close()
        h5file = None

    def create_next_run_group(self, h5file):
        return create_next_run_group(h5file, self.hdf5_base_loc)

    def save_results(self, hdf5_fn='adaptive_saw.h5',
                     res_num=None, force=False, other_specs = None):
        h5file = t.openFile(hdf5_fn, mode='a')

        if res_num is None:
            run_grp_loc, res_num = self.create_next_run_group(h5file)
        else:
            run_grp_loc = self.get_hdf5_location(res_num)

            try:
                run_grp = h5file.getNode(run_grp_loc)
                if force:
                    run_grp._f_remove(True)
                else:
                    raise Exception("Results already exist in " + run_grp)
            except t.NoSuchNodeError as e:
                pass # Node doesn't exist - don't need to delete
            
            h5file.createGroup(self.hdf5_base_loc, 'r' + res_num,
                               '', None, True)

        print "saving results for test in", run_grp_loc
        
        map(lambda name:
            h5file.createArray(run_grp_loc, name, getattr(self, name)),
            self.saved_res_fields)

	if other_specs is not None:
	    map(lambda key:
            h5file.createArray(run_grp_loc, key, other_specs[key]),
            other_specs.keys())

        self.res_num = res_num

        h5file.flush()
        h5file.close()
        h5file = None

    def get_num_results(self, hdf5_fn='adaptive_saw.h5'):
        h5file = t.openFile(hdf5_fn, mode='r')
        
        nresults = get_num_run_groups(h5file, self.hdf5_base_loc)
        
        h5file.flush()
        h5file.close()
        h5file = None

        return nresults

    def clear_all_results(self, hdf5_fn='adaptive_saw.h5'):
        h5file = t.openFile(hdf5_fn, mode='a')

        try:
            res_grp = h5file.getNode(self.hdf5_base_loc)
            res_grp._f_remove(True)
        except t.NoSuchNodeError as e:
            pass # Node does not exist - nothing to remove
        
        h5file.flush()
        h5file.close()
        h5file = None

    def plot_mh_ratio_moving_avg(self):
        plot_name = "mh_ratio_moving_avg"
        import utils
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(utils.compute_moving_average(self.MH_ratios))
        ax.set_ylim(0, 1)
        ax.set_xlim(0, len(self.MH_ratios))
        self.save_current_plot(plot_name)

    def plot_E_samples(self):
        plot_name = "E_samples"

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.E_samples)

        ax.set_xlabel('Steps')
        ax.set_ylabel('Energy')

        self.save_current_plot(plot_name)

    def plot_acf(self, max_lag = 500, start = 0):
        plot_name = "acf"

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(map(lambda i: utils.compute_acf(self.E_samples[start:], i),
                    xrange(max_lag+1)))
        ax.set_ylim(0, 1)
        ax.set_xlim(0, max_lag)
        self.save_current_plot(plot_name)

    def save_current_plot(self, plot_name):
        plot_dir = "plots/%s" % (self.graph.name)
        name = self.name
        if hasattr(self, 'nadaptations') and self.nadaptations is not None:
	    name = self.name + '_a' + str(self.nadaptations)
	if hasattr(self, 'nadaptations') and self.policy_size is not None:
	    name = name + '_p' + str(self.policy_size)
	
        plot_fn = "%s/%s_%s_r%d.pdf" % \
                  (plot_dir, plot_name, name, self.res_num)
        self.save_current_plot_to_file(plot_fn)

    def save_current_plot_to_file(self, plot_fn):
        import utils
        import os.path
        utils.ensure_dir_exists(os.path.dirname(plot_fn))

        import matplotlib.pyplot as plt
        plt.savefig(plot_fn, format='pdf')
        print "Saved plot to %s" % (plot_fn)

class SAW_Test(Sampler_Test):
    SAW_length_min = None
    SAW_length_max = None
    
    saved_res_fields = \
        set([]) | \
        Sampler_Test.saved_res_fields

class IM_Test(SAW_Test):
    S_ref = None
    
    saved_res_fields = \
        set(['SAW_length_max', 'SAW_length_min', 'SAW_lengths',
             'S_ref', 'used_gammas']) | \
        SAW_Test.saved_res_fields

    def __init__(self, graph, nmoves, SAW_length_min, SAW_length_max,
                 S_ref=None, S_in=None):
        im.load_im_sampler_lib()
        
        self.graph = graph
        self.nmoves = nmoves
        self.SAW_length_min = SAW_length_min
        self.SAW_length_max = SAW_length_max
        
        self.hdf5_base_loc = "/" + graph.name + "/" + self.name

        # Default reference state is ground state
        if S_ref is None:
            self.S_ref = self.graph.get_default_reference_state()
        else:
            self.S_ref = S_ref

        # Default initial state is reference state with half of its
        # bits flipped
        if S_in is None:
            self.S_in = self.graph.get_default_initial_state()
        else:
            self.S_in = S_in

class SARD_Test(SAW_Test):
    saved_res_fields = \
        set(['gamma_hi_bounds', 'gamma_lo_bounds',
             'iters_per_move_bounds',
             'max_SAW_length_bounds', 'min_SAW_length_bounds',
             'used_sigma_lengths', 'used_rho_lengths',
             'used_gammas_his', 'used_gammas_los',
             'used_iters_per_move', 'move_types',
             'used_P_LL', 'used_P_LH', 'used_P_HL']) | \
        SAW_Test.saved_res_fields

    def __init__(self, graph, nmoves,
                 gamma_hi_bounds, gamma_lo_bounds,
                 iters_per_move_bounds,
                 max_SAW_length_bounds, min_SAW_length_bounds,
                 strategy, P_LL = None, P_HL = None,
                 S_in=None, return_samples = False):
        sard.load_sard_sampler_lib()

	self.return_samples = return_samples

        self.graph = graph
        self.nmoves = nmoves
        
        self.gamma_hi_bounds = gamma_hi_bounds
        self.gamma_lo_bounds = gamma_lo_bounds
        self.max_SAW_length_bounds = max_SAW_length_bounds
        self.min_SAW_length_bounds = min_SAW_length_bounds
        self.iters_per_move_bounds = iters_per_move_bounds
        
        self.hdf5_base_loc = "/" + graph.name + "/" + self.name
        
        self.strategy = strategy
	
	self.P_HL = P_HL
	self.P_LL = P_LL
	
        # Default initial state is reference state with half of its
        # bits flipped
        if S_in is None:
            self.S_in = self.graph.get_default_initial_state()
        else:
            self.S_in = S_in

class IM_Unif_Test(IM_Test):
    name = "IM_Unif"

    def run(self):
        nruns_per_param_update = 1
        gammas = arange(0,2.0001,0.2) * self.graph.beta_true
        
        h_effective_in = self.graph.get_effective_fields(self.S_in)

        self.S_end, self.E_samples, self.E_proposeds, \
        self.log_f_fwds, self.log_f_revs, \
        self.MH_ratios, self.SAW_lengths, self.used_gammas \
        = im.IMRun_Unif\
          (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
           self.graph.edges_i, self.graph.J, self.graph.h,
           self.S_ref, self.S_in, h_effective_in,
           self.graph.beta_true, gammas,
           self.SAW_length_min, self.SAW_length_max,
           self.nmoves, nruns_per_param_update)

class SARD_Unif_Test(SARD_Test):
    name = "SARD_Unif"
    
    saved_res_fields = \
        set(['return_samples']) | \
        SARD_Test.saved_res_fields

    def run(self):
        nruns_per_param_update = 1
        h_effective_in = self.graph.get_effective_fields(self.S_in)
        
        
        if not self.return_samples:
	    self.S_end, self.E_samples, self.E_proposeds, \
	    self.log_f_fwds, self.log_f_revs, self.MH_ratios, \
	    self.move_types, \
	    self.used_sigma_lengths, self.used_rho_lengths, \
	    self.used_iters_per_move, \
	    self.used_gammas_his, self.used_gammas_los, \
	    self.used_P_LL, self.used_P_LH, self.used_P_HL \
	    = sard.SARDRun_Unif\
	    (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
	    self.graph.edges_i, self.graph.J, self.graph.h,
	    self.S_in, h_effective_in, self.graph.beta_true,
	    self.gamma_hi_bounds, self.gamma_lo_bounds,
	    self.iters_per_move_bounds,
	    self.max_SAW_length_bounds, self.min_SAW_length_bounds,
	    self.nmoves, nruns_per_param_update, self.return_samples, self.strategy, self.P_LL, self.P_HL)
	else:
	    return sard.SARDRun_Unif\
	    (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
	    self.graph.edges_i, self.graph.J, self.graph.h,
	    self.S_in, h_effective_in, self.graph.beta_true,
	    self.gamma_hi_bounds, self.gamma_lo_bounds,
	    self.iters_per_move_bounds,
	    self.max_SAW_length_bounds, self.min_SAW_length_bounds,
	    self.nmoves, nruns_per_param_update, self.return_samples, self.strategy, self.P_LL, self.P_HL)

class Kawasaki_Sampler_Test(IM_Test):
    name = "Kawasaki_Sampler"

    def __init__(self, graph, nmoves, S_ref=None, S_in=None):
        super(Kawasaki_Sampler_Test, self).__init__\
            (graph, nmoves, 0, 0, S_ref, S_in)

    def run(self):
        im.load_im_sampler_lib()
        
        h_effective_in = self.graph.get_effective_fields(self.S_in)

        self.S_end, self.E_samples, self.E_proposeds, \
        self.log_f_fwds, self.log_f_revs, \
        self.MH_ratios, self.SAW_lengths, self.used_gammas \
        = im.Kawasaki_Sampler\
          (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
           self.graph.edges_i, self.graph.J, self.graph.h,
           self.S_ref, self.S_in, h_effective_in,
           self.graph.beta_true, self.nmoves)

class IM_Expert_Test(IM_Test):
    name = "IM_Expert"

    def run(self):
        nruns_per_param_update = 1
        gamma = self.graph.beta_true # gamma is set to beta_true

        h_effective_in = self.graph.get_effective_fields(self.S_in)

        self.S_end, self.E_samples, self.E_proposeds, \
        self.log_f_fwds, self.log_f_revs, \
        self.MH_ratios, self.SAW_lengths, self.used_gammas \
        = im.IMRun_Expert\
          (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
           self.graph.edges_i, self.graph.J, self.graph.h,
           self.S_ref, self.S_in, h_effective_in,
           self.graph.beta_true, gamma,
           self.SAW_length_min, self.SAW_length_max,
           self.nmoves, nruns_per_param_update)

class IM_EGO_Test(IM_Test):
    name = "IM_EGO"

    nadapted_moves = None
    nruns_per_param_update = None
    nadaptations = None

    gamma_min = None
    gamma_max = None

    saved_res_fields = \
        set(['nadapted_moves', 'nruns_per_param_update',
             'nadaptations', 'gp_X', 'gp_Y', 'rewards',
             'kernel_hyperparms', 'gp_noise',
             'gamma_min', 'gamma_max', 'policy']) | \
        IM_Test.saved_res_fields

    def __init__(self, graph,
                 nmoves, nadaptations, nruns_per_param_update,
                 policy_size,
                 SAW_length_min, SAW_length_max, strategy, exp_temp, 
                 MH_rate, S_ref=None, S_in=None):
        super(IM_EGO_Test, self).__init__\
            (graph, nmoves,
             SAW_length_min, SAW_length_max, S_ref, S_in)

	self.exp_temp = exp_temp
	self.MH_rate = MH_rate
	self.strategy = strategy
        self.policy_size = policy_size
        self.hdf5_base_loc += "/p" + str(self.policy_size)

        self.nadapted_moves = \
            array([get_num_adapted_moves(nruns_per_param_update,
                                         nadaptations)])
        self.nadaptations = nadaptations
        self.nruns_per_param_update = nruns_per_param_update

        self.gammas = arange(0,2.0001,0.2) * self.graph.beta_true
        self.h_effective_in = self.graph.get_effective_fields(self.S_in)

        self.sampler_fn = \
            lambda : im.IMRun_EGO\
            (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
             self.graph.edges_i, self.graph.J, self.graph.h,
             self.S_ref, self.S_in, self.h_effective_in,
             self.graph.beta_true, self.gammas,
             self.SAW_length_min, self.SAW_length_max,
             self.nmoves, self.nadaptations,
             self.nruns_per_param_update, self.policy_size,
             self.strategy, self.exp_temp, self.MH_rate)

    def run(self):
        self.S_end, self.E_samples, self.E_proposeds, \
        self.log_f_fwds, self.log_f_revs, \
        self.MH_ratios, self.SAW_lengths, self.used_gammas, \
        self.rewards, gp, self.gamma_min, self.gamma_max, \
        self.policy \
        = self.sampler_fn()

        self.gp_X = gp.X
        self.gp_Y = gp.Y

        self.kernel_hyperparms = gp.kernel.getHyperparams()
        self.gp_noise = gp.noise

    def load_results(self, *args):
        super(IM_EGO_Test, self).load_results(*args)

        from ego.gaussianprocess import GaussianProcess
        from ego.gaussianprocess.kernel import GaussianKernel_ard

        # Build the gp
        kernel = GaussianKernel_ard(self.kernel_hyperparms)
        self.gp = GaussianProcess(kernel, noise=self.gp_noise)

        for x, y in zip(self.gp_X, self.gp_Y):
            self.gp.addData(x, y)

    def save_current_plot(self, plot_name):
        plot_dir = "plots/%s" % (self.graph.name)
        plot_fn = "%s/%s_%s_p%d_r%d.pdf" % \
                  (plot_dir, plot_name, self.name,
                   self.policy_size, self.res_num)
        self.save_current_plot_to_file(plot_fn)


    def compute_gp_posterior_over_grid(self):
        x1 = arange(self.SAW_length_min, self.SAW_length_max, 1, dtype=float)
        x2 = arange(self.gamma_min, self.gamma_max,
                    (self.gamma_max - self.gamma_min) / 100,
                    dtype=float)
        
        X1, X2 = meshgrid(x1, x2)

        m = zeros_like(X1)
        v = zeros_like(X1)

        for i in xrange(X1.shape[0]):
            for j in xrange(X1.shape[1]):
                z = array([X1[i, j], X2[i, j]])

                res = self.gp.posterior(z)
                m[i, j] = res[0]
                v[i, j] = res[1]

        return X1, X2, m, v

    def plot_surface_on_grid(self, X1, X2, s):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        
        fig = plt.figure()
        ax = Axes3D(fig, elev=70)
        surf = ax.plot_surface(X1, X2, s, rstride=1, cstride=1, cmap=cm.jet,
                       linewidth=0, antialiased=True)
        ax.set_zlim3d(-1.01, 1.01)

        ax.w_zaxis.set_major_locator(LinearLocator(10))
        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.03f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

    def plot_gp_mean(self):
        plot_name = "gp_mean"
        
        X1, X2, m, _ = self.compute_gp_posterior_over_grid()
        self.plot_surface_on_grid(X1, X2, m)
        self.save_current_plot(plot_name)

    def plot_flat_gp_mean(self):
        plot_name = "flat_gp_mean"

        import matplotlib.pyplot as plt
        from matplotlib import cm

        _, _, m, _ = self.compute_gp_posterior_over_grid()
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(m, cmap=cm.jet, interpolation='nearest', origin='lower',
                  extent=[self.SAW_length_min, self.SAW_length_max,
                          self.gamma_min, self.gamma_max],
                  aspect=self.SAW_length_max/self.gamma_max)
        ax.set_xlabel('SAW length $k$')
        ax.set_ylabel('Energy-biasing parameter $\gamma$')
        
        self.save_current_plot(plot_name)

    def plot_average_gp_mean(self):
        plot_name = "average_gp_mean"

        nresults = self.get_num_results()
        
        X1 = None
        X2 = None
        m_avg = None

        if self.res_num is not None:
            old_res_num = self.res_num
        else:
            old_res_num = 0

        for i in xrange(nresults):
            self.load_results(i)
            X1, X2, m, _ = self.compute_gp_posterior_over_grid()

            if m_avg is None:
                m_avg = m
            else:
                m_avg += m

        m_avg /= nresults

        self.plot_surface_on_grid(X1, X2, m_avg)
        self.save_current_plot(plot_name)

        self.load_results(old_res_num)

    def plot_flat_average_gp_mean(self):
        plot_name = "flat_average_gp_mean"

        import matplotlib.pyplot as plt
        from matplotlib import cm

        nresults = self.get_num_results()

        m_avg = None

        if self.res_num is not None:
            old_res_num = self.res_num
        else:
            old_res_num = 0
            
        for i in xrange(nresults):
            self.load_results(i)
            _, _, m, _ = self.compute_gp_posterior_over_grid()

            if m_avg is None:
                m_avg = m
            else:
                m_avg += m

        m_avg /= nresults

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(m_avg, cmap=cm.jet, interpolation='nearest', origin='lower',
                  extent=[self.SAW_length_min, self.SAW_length_max,
                          self.gamma_min, self.gamma_max],
                  aspect=self.SAW_length_max/self.gamma_max)
        ax.set_xlabel('SAW length $k$')
        ax.set_ylabel('Energy-biasing parameter $\gamma$')

        self.save_current_plot(plot_name)

        self.load_results(old_res_num)

    def plot_gp_var(self):
        plot_name = "gp_var"
        
        X1, X2, _, v = self.compute_gp_posterior_over_grid()
        self.plot_surface_on_grid(X1, X2, v)
        self.save_current_plot(plot_name)

    def plot_gp_query_points(self):
        plot_name = "gp_query_points"
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(self.gp.X[:, 0], self.gp.X[:, 1], c=self.gp.Y)

        self.save_current_plot(plot_name)

    def plot_parameter_cdf(self):
        pass

    def make_policy_dict(self):
        policy = map(lambda x, y: (round(x[0],3), round(y[0], 3)),
                     self.policy[::2], self.policy[1::2])
        policy_dict = {}
        for x in policy:
            policy_dict[x] = policy_dict.get(x, 0) + 1

        return policy_dict

    def plot_policy_as_histogram(self):
        plot_name = "policy_histogram"

        import matplotlib.pyplot as plt
        
        xpos = []
        ypos = []
        dz = []
        
        policy_dict = self.make_policy_dict()
        for k, v in policy_dict.iteritems():
            xpos.append(k[0])
            ypos.append(k[1])
            dz.append(v)

        dy_size = (self.gammas[-1] - self.gammas[0]) / 100

        zpos = zeros_like(xpos)
        dx = ones_like(xpos)
        dy = dy_size * ones_like(ypos)

        colour_max = max(dz)
        colours = map(lambda x: plt.cm.jet(x / colour_max), dz)

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        
        fig = plt.figure()
        ax = Axes3D(fig, elev=70)

        ax.w_zaxis.set_major_locator(LinearLocator(10))
        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.03f'))

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,
                 edgecolors=(0.0, 0.0, 0.0, 0.0),
                 color=colours)
        ax.set_xlim3d(self.SAW_length_min, self.SAW_length_max)
        ax.set_ylim3d(self.gammas[0], self.gammas[-1])

        self.save_current_plot(plot_name)

    def plot_policy_as_flat_surface():
        plot_name = 'flat_policy'

class IM_EGO_premixed_Test(IM_EGO_Test):
    name = "IM_EGO_premixed"

    def __init__(self, graph,
                 nmoves, nadaptations, nruns_per_param_update,
                 policy_size,
                 SAW_length_min, SAW_length_max,
                 S_ref=None, S_in=None):
        super(IM_EGO_premixed_Test, self).__init__\
            (graph, nmoves, nadaptations, nruns_per_param_update,
             policy_size,
             SAW_length_min, SAW_length_max, S_ref, S_in)

        self.sampler_fn = \
            lambda : im.IMRun_EGO_premixed\
            (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
             self.graph.edges_i, self.graph.J, self.graph.h,
             self.S_ref, self.S_in, self.h_effective_in,
             self.graph.beta_true, self.gammas,
             self.SAW_length_min, self.SAW_length_max,
             self.nmoves, self.nadaptations,
             self.nruns_per_param_update, self.policy_size)

class SARD_EGO_Test(SARD_Test):
    name = "SARD_EGO"

    nadapted_moves = None
    nruns_per_param_update = None
    nadaptations = None

    saved_res_fields = \
        set(['nadapted_moves', 'nruns_per_param_update',
             'nadaptations', 'gp_X', 'gp_Y', 'rewards',
             'kernel_hyperparms', 'gp_noise', 'policy']) | \
        SARD_Test.saved_res_fields

    def __init__(self, graph,
                 nmoves, nadaptations, nruns_per_param_update,
                 policy_size,
                 gamma_hi_bounds, gamma_lo_bounds,
                 iters_per_move_bounds, strategy, exp_temp, MH_rate,
                 max_SAW_length_bounds, min_SAW_length_bounds,
                 S_in=None, return_samples=1):
        super(SARD_EGO_Test, self).__init__\
            (graph, nmoves, gamma_hi_bounds, gamma_lo_bounds,
             iters_per_move_bounds,
             max_SAW_length_bounds, min_SAW_length_bounds, S_in)

	self.max_SAW_length_bounds = [self.max_SAW_length_bounds[i] - \
	    self.min_SAW_length_bounds[i] for i in range(len(self.min_SAW_length_bounds))]

	self.gamma_hi_bounds = [self.gamma_hi_bounds[i] - \
	    self.gamma_lo_bounds[i] for i in range(len(self.gamma_lo_bounds))]

        self.policy_size = policy_size
        self.hdf5_base_loc += "/p" + str(self.policy_size)

        self.nadapted_moves = \
            array([get_num_adapted_moves(nruns_per_param_update,
                                         nadaptations)])
        self.nadaptations = nadaptations
        self.nruns_per_param_update = nruns_per_param_update

        h_effective_in = self.graph.get_effective_fields(self.S_in)

        self.sampler_fn = \
            lambda : sard.SARDRun_EGO \
            (self.graph.nbrs_list_i, self.graph.incident_edge_list_i,
             self.graph.edges_i, self.graph.J, self.graph.h,
             self.S_in, h_effective_in, self.graph.beta_true,
             self.gamma_hi_bounds, self.gamma_lo_bounds,
             self.iters_per_move_bounds,
             self.max_SAW_length_bounds, self.min_SAW_length_bounds,
             self.nmoves, self.nadaptations,
             self.nruns_per_param_update, self.policy_size, strategy, exp_temp, MH_rate, return_samples)
    def make_policy_dict(self):
	reshaped = self.policy.reshape(8, self.policy_size)
	
        policy = map(lambda x, y: (round(x[0],0), round(y[0], 0)),
                     self.policy[1::8], self.policy[1::8]+self.policy[2::8])
                     
        policy_dict = {}
        for x in policy:
            policy_dict[x] = policy_dict.get(x, 0) + 1

        return policy_dict

    def plot_policy_as_histogram(self):
        plot_name = "policy_histogram"

        import matplotlib.pyplot as plt
        
        xpos = []
        ypos = []
        dz = []
        
        policy_dict = self.make_policy_dict()
        for k, v in policy_dict.iteritems():
            xpos.append(k[0])
            ypos.append(k[1])
            dz.append(v)

	dx_size = (self.min_SAW_length_bounds[1] - self.min_SAW_length_bounds[0])/200
	dy_size = (self.max_SAW_length_bounds[1] - self.max_SAW_length_bounds[0])/200 + dx_size
	
	
        zpos = zeros_like(xpos)
        dx = dx_size * ones_like(xpos)
        dy = dy_size * ones_like(ypos)

        colour_max = max(dz)
        colours = map(lambda x: plt.cm.jet(x / colour_max), dz)

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        
        fig = plt.figure()
        ax = Axes3D(fig, elev=70)

        ax.w_zaxis.set_major_locator(LinearLocator(10))
        ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.03f'))

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz,
                 edgecolors=(0.0, 0.0, 0.0, 0.0),
                 color=colours)
	ax.set_xlim3d(self.min_SAW_length_bounds[0], self.min_SAW_length_bounds[1])
	ax.set_ylim3d(self.min_SAW_length_bounds[0] + self.max_SAW_length_bounds[0], self.min_SAW_length_bounds[1] + self.max_SAW_length_bounds[1])

        self.save_current_plot(plot_name)

    def run(self):
        self.S_end, self.E_samples, self.E_proposeds, \
        self.log_f_fwds, self.log_f_revs, \
        self.MH_ratios, self.move_types, \
        self.used_sigma_lengths, self.used_rho_lengths, \
        self.used_iters_per_move, \
        self.used_gammas_his, self.used_gammas_los, \
        self.used_P_LL, self.used_P_LH, self.used_P_HL, \
        self.rewards, gp, self.policy, self.samples \
        = self.sampler_fn()

        self.gp_X = gp.X
        self.gp_Y = gp.Y

        self.kernel_hyperparms = gp.kernel.getHyperparams()
        self.gp_noise = gp.noise
    
    def load_results(self, *args):
        super(SARD_EGO_Test, self).load_results(*args)

        from ego.gaussianprocess import GaussianProcess
        from ego.gaussianprocess.kernel import GaussianKernel_ard

        # Build the gp
        kernel = GaussianKernel_ard(self.kernel_hyperparms)
        self.gp = GaussianProcess(kernel, noise=self.gp_noise)

        for x, y in zip(self.gp_X, self.gp_Y):
            self.gp.addData(x, y)

	print shape(self.S_end)
	print shape(self.E_samples)
	print shape(self.MH_ratios)
	print shape(self.E_proposeds)
