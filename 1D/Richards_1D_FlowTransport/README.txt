This folder contains the following Matlab scripts and files:


'main_Richiards_1D_FlowTransp_benchmark.m' is the code which solves the 1D coupled flow and 
	transport problem for the 1D setup inspired by the 2D benchmark problem considered 
	by [Schneid, 2000; List and Radu, 2016].

'BGRW_1D_Alt_rand_benchmark.m' is the function solving the transport step of the coupled 
	problem by the BGRW algorithm.

'K_r.m' is a function providing lognormally distributed random hydraulic conductivity fields.

'Kraichnan_Exp_param.m' and 'Kraichnan_Gauss_param.m' are functions used to compute random 
	phases and wave numbers used to generate lognormal random fields.

'initstate.m' is used to fix the seed of the random number generator.

'theta_GM.m' is the function for the calculation of the unsaturated/saturated water content 
	according to van Genuchten-Mualem parameterization.

'plot_fig.m' is a function providing pictures with results of the computation.

'comparison_clay_loam_1D.m' is a code comparing the numerical solutions 'clay'/'clay_det' and 
	'loam'/'loam_det' of the benchmark problem for random/determinstic sauturated-hydraulic 
	conductivity.

'convf.mat' and 'convt.mat' are temporary files containing norms of the convergence criterium for 
	the flow and transport solvers, respectively.