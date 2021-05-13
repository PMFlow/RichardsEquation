## Data repository accompanying the paper "Global random walk solvers for fully coupled flow and transport in saturated/unsaturated porous media" by Nicolae Suciu, Davide Illiano, Alexander Prechtel, Florin A. Radu
This is a data repository accompanying the paper "Global random walk solvers for fully coupled flow and transport in saturated/unsaturated porous media" by Nicolae Suciu, Davide Illiano, Alexander Prechtel, Florin A. Radu (https://doi.org/10.1016/j.advwatres.2021.103935)

The repository contains Matlab codes and benchmark tests for the Global Random Walk (GRW) approach to solve coupled nonlinear problems of flow and transport in porous media.

The GRW solutions for a variety of one- and two-dimensional problems are grouped together in the folders 1D and 2D as follows.
#
### 1D/Richards_1D
- illustrates the construction of the GRW solution, the L-scheme used to linearize the nonlinear Richards equation, and the transition from unsaturated to saturated flow regime, for exponential parameterization of the water content and hydraulic conductivity as a function of pressure head.
### 1D/Richards_1D_code_verification
- accuracy and order of convergence assessments using analytical manufactured solutions.
### 1D/Richards_1D_validation_tests
- comparison with experimental data and exact analytical solutions.
### 1D/Richards_1D_FlowTransport
- benchmark for coupled flow and transport in unsaturated/saturated soils with van Genuchtem-Mualem parameterization and random saturated-conductivity.
##
### 2D/Richards_2D
- flow benchmark for clay- and loam-soil models with van Genuchtem-Mualem parameterization and random saturated-conductivity.
### 2D/Richards_2D_FlowTransport
- benchmark for coupled flow and transport in unsaturated/saturated soils with van Genuchtem-Mualem parameterization and random saturated-conductivity.
### 2D/Saturated_2D_FlowTransport
- flow in saturated aquifers at regional scale with random hydraulic conductivity, flow with random recharge of the aquifer system, and Monte Carlo simulations advection-dispersion in random velocity fields.
