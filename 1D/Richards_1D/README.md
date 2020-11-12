This folder contains the following files:
#

- 'main_L_GRW_Richy_1D.m' is the Matlab code for the L_GRW scheme used to solve flow-problems for
	two numerical scenarios considered in 'Richy-1D' software to obtain reference solutions.
	
- 'theta_exp.m' provides the unsaturated/saturated water content as a function of pressure head 
	according to the exponential parameterization.
	
- 'IC_L_GRW_Richy_1D.m' computes the initial conditions 'IC_GRW_...' as solutions of steady-state 
	problems designed for the two scenarios.
	
- 'Richy_...' are files containing the Richy-1D reference solutions and 'IC_Richy_...' are files 
	containing the corresponding initial conditions.
	
- 'solution_scenario1.mat' and 'solution_scenario2.mat' contain the numerical L_GRW solutions. 

- 'dt_scenario1' and 'dt_scenario2' are files which store the variable time step updated by the 
	iterative procedure used in the 'L_GRW' scheme.
	
- 'comparison_IC_GRW_Richy.mm' and 'comparison_solution_GRW_Richy.m' are codes used to compare  the 		
	initial conditions and the solutions.
