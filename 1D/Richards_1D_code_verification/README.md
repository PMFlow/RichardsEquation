This folder contains the following Matlab scrips and files:
#

- 'main_Richy1D_FlowTransp_test.m' is the main porgram which solves the coupled flow and 
	transport problem and compares the numerical solution to the analytical reference 
	solution p(t,x)=-tx(1-x)-1, c(t,x)=tx(1-x)+1.

- 'BGRW_1D_Alt_rand_benchmark.m' is the function solving the transport step of the coupled 
	problem by the BGRW algorithm.

- 'theta.m' is the function for water content in the degenerate problem.

- 'Fsg.m' is the function for the source terms in the degenerate eproblem.

- 'plot_conv_Alt.m' is a function which plotes the norm of the convergence criterium versus 
	the number of iterations at different times moments.

- 'plot_fig.m' is a function plotting the numerical solution at the final time.

- 'test_dt' is a file containing the time step as a function of the current computation time.

- 'convf.mat' and 'convf.mat' are fiels containing norms of the convergence criterium for the 
	flow and transport solves.

