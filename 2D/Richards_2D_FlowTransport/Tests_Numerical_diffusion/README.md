This folder contains the following Matlab scripts:
#

- 'main_tests_num_diff.m' is the main porgram which computes the constant velocity and calls 
	either 'testBGRW_const.m' or 'testGRW_const.m' to estimate numerical diffusion errors.

- 'testBGRW_const.m' is a function which solves the transport problem with constant coefficients by a 
	Biased_GRW algorithm and estimates the numerical diffusion by comparing numerically estimated 
	diffusion coefficients to the nominal coefficients of the analytical Gaussian distribution.

- 'testGRW_const.m' is a function performing the same operations as above with an unbiased GRW algorithm.

- 'Gauss_IC.m' is the function providing the theoretical Gaussian solution of the advection-diffusion 
	equation with constant coefficients.
