This folder contains the following files:
#

- 'main_Richards1D_Warrick.m' is the Matlab code for the GRW scheme to solve the 1D Richards 
	equation in mixed-formulation with numerical parameters from [Warrick et al. 1985].

- 'theta_GM.m' is the function for the calculation of the water content \theta according to 
	the van Genuchten-Mualem parameterization.

- 'test_dt' is a file containing the variable time steps.

- 'convf' is a file containing norms of the convergence criterium at times t<=T.

- 'plot_conv.m' is the function plotting the convergence norms.

- 'pT', 'thT' and 'qT' are files containing solutions at different final times T.

- 'plot_T.m" is the function plotting solutions at the final times T.

- 'Warricks_solution_n1_5.m' is a Matlab code to compute Warrick's solutions for n=1.5 at different 
	times T according to [Warrick et al., 1985, Table 3].

- 'WX_n1_5.mat' is the file containing the Warrick's solution.

- 'errors_X.m' is the Matlab code which interpolates the GRW solutions \theta(z) from 'thT' and 
	computes relative errors with respect to Warrick's solution 'WX_n1_5.mat'.
