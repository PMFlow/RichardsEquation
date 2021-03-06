This folder contains the following Matlab scripts and files:
#

- 'main_ensemble_coefficients.m' is a 2D  GRW algorithm for advection-diffusion in random fields using velocity realizations stored in the sub-folder 'GAUSS'; estimated mean velocities and ensemble 		dispersion coefficients are stored in 'ensemble_coefficients.mat'.

- 'main_realizations.m' computes random velocity realizations stored in the sub-folder 'GAUSS' by utilizing the flow solver 'realiz_Gauss_GRW.m'.

- 'convergence_test_GRW_Flow.m' analyzes the convergence of the flow solver.

- 'Kraichnan_Gauss_param.m', 'K_repmat.m', 'Gauss_IC.m' are functions used in the 'main' codes.

- 'linear approximation' folder contains codes to compute the 1-st order (linear) approximations of	ensemble dispersion coefficients stored in 'linear_approximation.mat'

- 'comparison' compares GRW coefficients stored in 'ensemble_coefficients.mat' to reference 1-st order coefficients from 'linear_approximation.mat'.
