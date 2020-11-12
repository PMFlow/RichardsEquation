This folder contains the following Matlab scripts and files:
#

- 'main_Richiards_2D_FlowTransp_benchmark.m' it the code which solves the 2D problem obtained by coupling a nonlinear reactive-transport problem with the 2D benchmark flow-problem considered by [Schneid, 2000; List and Radu, 2016].

- 'BGRW_Alt_rand_benchmark.m' is the function solving the transport step of the coupled problem by the BGRW algorithm.

- 'K_r.m' is a function providing lognormally distributed random hydraulic conductivity fields.

- 'Kraichnan_Exp_param.m' and 'Kraichnan_Gauss_param.m' are function used to compute random 
	phases and wave numbers used to generate lognormal random fields.

- 'initstate.m' is used to fix the seed of the random number generator.

- 'theta_GM.m' is the function for the calculation of the unsaturated/saturated water content 
	according to van Genuchten-Mualem parameterization.

- 'plot_fig.m' is a function providing three-dimensional plots of the solution at the 
	final time.

- 'convergence_plots.m' is a code providing plots of the convergence norms for the flow and 
	transport solvers.

- 'contour_plots.m' is a code providing plane contours of the solution at the final time.

- 'convf_...' and 'convt_...' are mat-files containing norms of the convergence criterium for 
	loam- and clay-solutions.

- 'solution_loam.mat' and 'solutions_clay.mat are the files containing the solution of the 
	coupled flow and transport problems for loam and clay soil models, respectively.
