This folder contains the following Matlab scripts and files:
#

- 'main_Richiards_2D_Flow_benchmark.m' is the code which solves the 2D benchmark flow-problem 
	considered by [Schneid, 2000; List and Radu, 2016].

- 'Richy_2D_Flow_benchmark_IC_SongJiang.m' is the code which reproduces the initial conditions
  	fabricated by Song and Jiang (2023, Figs. 1 and 2, https://checlams.github.io/assets/pdf/paper/songescape2023.pdf)

- 'theta_GM.m' is the function for the calculation of the unsaturated/saturated water content 
	according to van Genuchten-Mualem parameterization.

- 'plot_fig.m' and 'plot_contours.m' are functions providing three-dimensional plots and plane 
	contours of the solution at the final time.

- 'sol_GRW_2D_loam.mat' and 'sol_GRW_2D_clay.mat' are files containing the solutions for laom 
	and clay soil models, respectively.
