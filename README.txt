This package includes code for the SARDONICS sampler, the Intracluster Moves sampler as well as code for adapting these samplers by using Bayesian Optimization (aka GP-bandits).


===================================
BACKGROUND:
===================================
For detailed background information, please refer to the 2 references cited below. This package also contains the code of a modified version of the Swendsen-Wang sampler by Iain Murray. For detail please refer to the folder 'swendsen_wang'.


===================================
INSTALLATION:
===================================

This package was compiled on a Linux machine. But it should be easily modifiable to run on a Macintosh. This package requires several python packages including 'numpy', 'scipy', 'matplotlib', and 'pytalbes'. For the ease of installation, it is recommanded to install enthought which includes all the aformentioned packages. For simple installation, just run setup.sh. In case of an compilation error, please go through setup.sh and modify the relevant make files as needed. The code contains both c and python. If there is any problem compiling the code, please contact the author Ziyu Wang at ziyuw@cs.ubc.ca


===================================
KNOW ROUGH EDGES:
===================================
Currently, due to memory allocation issues, the code crashes on some graphs sizes. 



===================================
REFERENCES:
===================================

@article{hamze2011self,
  title={Self-Avoiding Random Dynamics on Integer Complex Systems},
  author={Hamze, F. and Wang, Z. and de Freitas, N.},
  journal={Arxiv preprint arXiv:1111.5379},
  year={2011}
}

@article{mahendran2011bayesian,
  title={Bayesian optimization for adaptive MCMC},
  author={Mahendran, N. and Wang, Z. and Hamze, F. and de Freitas, N.},
  journal={Arxiv preprint arXiv:1110.6497},
  year={2011}
}