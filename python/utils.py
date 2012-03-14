# default_mpl_backend = None # specified by .matplotlibrc
default_mpl_backend = 'Agg'

def load_matplotlib_backend(backend=default_mpl_backend):
    if backend is None:
        return
    
    import matplotlib
    #matplotlib.use(backend)
    
import numpy as np
def compute_moving_average(seq, mv_avg_window_size=100):
    # Compute the moving average of a sequence
    mv_avg_f = np.ones(mv_avg_window_size) / mv_avg_window_size
    avged_seq = np.convolve(mv_avg_f, np.resize(seq, (seq.shape[0],)))
    return avged_seq

def compute_acf(seq, k):
    # seq needs to be a numpy array, not a numpy matrix!
    n = len(seq)
    mu = np.mean(seq)
    sig_sqrd = np.var(seq)

    return sum((seq[:(n-k)] - mu) * (seq[k:n] - mu)) / ((n-k) * sig_sqrd)

class PathNotDirError(Exception):
    def __init__(self, path):
        self.path = path
    def __str__(self):
        return repr(self.path)

import os
def ensure_dir_exists(path):
    try:
        stat_res = os.stat(path)
        
        import stat
        mode = stat_res[stat.ST_MODE]
        
        if stat.S_ISDIR(mode):
            return
        else:
            raise PathNotDirError(path)
    except os.error as e:
        pass
    os.makedirs(path)



