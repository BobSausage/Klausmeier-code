# A_var

This folder contains all files relating to NSKM and SKM under changing rainfall

## solve_ivp_dkm_2D_A_var
This is the main file for the spatial model it contains the solve class. This class solves the system of ODEs in 2 dimensions.
An example of how to use the class:
```python
>>> b = solve(50,21, tmax = 200)      # this corresponds to a 50x50 matrix and will solve up to t=200 with 21 time steps including 0 and 200
>>> b.z0 = np.concatenate((read_ic('presets\\u8.txt'),np.ones(2500)*0.2)) # this gives some interesting initial conditions 
>>> b.solve2d(1,'BDF')             # this solves the system using BDF method with rainfall fn 1
effective A=4.300000612812226	   # it returns information on the average rainfall value and time taken to solve
time taken = 67.27607941627502
>>> b.visualise(3)  		   # returns a figure of the system at time interval 3 (in this case t=30)
```

Some parameters have been made class variables so they can be edited without having to re-initialise the class.
The important parameters are:
save : default set to 0, when set to 1 this allows to save plots to the plots folder
speed: default set to 0, when set to 1 this does not show diagrams ( making it easier to generate and save plots without having to stay at the computer)
amin : set by argument in initialising, parameter of rainfall
ka   : default set to 1, parameter of rainfall model 1&2
kb   : default set to 0.5, parameter of rainfall model 2
q    : default set to 1.3518, parameter of rainfall 1&2 (this should give A_eff = 1.1776 for rainfall 1)
r    : default set to 1, parameters of rainfall 1*2 (this should give A_eff = 1.1776 for rainfall 1)
z0  : default empty must be set to solve, this is the initial condition must have length 2*n^2 

This was used to generate Figures 6.4, 6.5, 6.6

## solve_ivp_dkm_2D_A_var_sloped
This contains the additional parameters and calculations to solve for sloped terrain or planes with water dispersal. It can be used for flat planes with no dispersal by setting parameters d and e to 0. Input is similar to above. 

This was used to generate Figures 6.7, 6.8, 6.9,6.10

## A_models
This file contains the rainfall functions. used for Figure 6.1

## jacobian_maker
This file contains the function intended to be used to generate the Jacobian for stiff methods. This was instead implemented directly into the main class file. 

## nskm_solve_ivp_varying_A
This file contains basic_fig_plot which was used for Figure 6.3 which shows the solution for given initial conditions over a time period. 
It also contains a function to find where in the rainfall cycle certain initial conditions would lead to late time plant density. 

## optimal_L
This file contains the optimised laplacian operators. 

## Final_phase_plane_special_plot_var_a
This plots the contours and basins of attraction. To plot onto already existing basin png the first lines of the from_wells function should be uncommented. The original code for the arrows and contours can be found in the "Other" folder labelled "online_python_phaseplane_plotter" however it was heavily edited by me for this file.  