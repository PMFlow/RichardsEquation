This folder contains the following files:


'main_Richiards_1D_drainage.m' is the Matlab code for the GRW scheme to solve the 1D Richards 
    equation for free drainage conditions (Hydrus 1D, 'drainage in a large caisson' example).

'theta_GM.m' is the function for the calculation of the water content \theta according to 
	the van Genuchten-Mualem parameterization.

'test_dt' is a file containing variable the time steps.

'convf' is a file containing norms of the convergence criterium at times t<=T.

'plot_conv.m' is the function plotting the convergence norms.

'pT', 'thT' and 'qT' are files containing GRW solutions at different final times T.

'plot_T.m" is the function plotting GRW solutions at the final times T.

'thH' and 'pH' are files containing Hydrus 1D solutions. 

'plot_theta.m' is the Matlab code for comparison of Experiment-Hydrus 1D-GRW solutions.