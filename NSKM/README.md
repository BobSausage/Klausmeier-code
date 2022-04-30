# NSKM
This folder contains code for the non-spatial model. (chapter 3)

## accuracy
This is plots viable values of A against B using mpmath to 1000dp accuracy. This was used to make Figure 3.1. It also contains the function getpoints which returns the saddle point and stable point for plant survival, given A and B.

## odeint_nskm
This contains code to solve the non-spatial model for given initial conditions using basic_fig_plot. It also contains functions to check many areas at once 
'roab' runs over values of a and b returns a plot for each pairing of and b values. 
'roz' runs over values of u and w returns a plot of the basins of attraction. This was used to generate Figure 3.3.

## Final_phase_plane_special_plot
This contains code for Figure 3.2. This plots the contours on basins of attraction. 
This contains two functions special_plot() which plots the contours on the basisn of attraction and from_pot() which plots the contours on an already existing png.
