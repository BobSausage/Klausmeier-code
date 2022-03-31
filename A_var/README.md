# A_var

This folder contains all files relating to NSKM and SKM under changing rainfall

## solve_ivp_dkm_2D_A_var
This is the main file for the spatial model it contains the solve class. This class solves the system of ODEs in 2 dimensions.
An example of how to use the class:
```python
>>> b = solve(50, tmax = 200)      # this corresponds to a 50x50 matrix and will solve over 200 years
>>> b.z0 = np.concatenate((u8,w2)) # this gives some interesting initial conditions 
>>> b.solve2d('BDF')               # this solves the system  
effective A=4.300000612812226	   # it returns information on the average rainfall value and time taken to solve
time taken = 67.27607941627502
>>> b.visualise(3)  		   # returns a figure of the system at time interval 3 (in this case 60 years)
```
At the bottom of the file are some initial conditions u1 - 13 and w1 - 2 to make writing initial conditions easier to set.

Some parameters have been made class variables so they can be edited without having to re-initialise the class.
The important parameters are:
save : default set to 0, when set to 1 this allows to save plots to the plots folder
speed: default set to 0, when set to 1 this does not show diagrams ( making it easier to generate and save plots without having to stay at the computer)
amin : set by argument in initialising, parameter of rainfall
kk   : default set to 2, parameter of rainfall
pq   : default set to 2, parameter of rainfall

## attempting_2D_matrix
This file contains the function which is called to generate the dispersal operator, it generates an n^2 x n^2 matrix with periodic boundary conditions.

## jacobian_maker
This file contains the function intended to be used to generate the Jacobian for stiff methods. This was instead implemented directly into the main class file. 
