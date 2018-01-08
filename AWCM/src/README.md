---
title: "README"
---

## Introduction
The following programs within the AWCM directory are numerical algorithms for solving partial 
differential equations. Specifically, the method used is the Adaptive Wavelet Collocation Method. This 
work was contributed to by the efforts cited in the References section. 

## File structure
Currently, the `src` directory contains a number of sub-directories and files necessary for compilation.

## Input file
The input file specifies all parameters which are input to the program. These include grid and refinement
parameters, physical parameters, numerical scheme choices, as well as options for plotting, producing output
and more. An example input file is shown below 

```{r,eval=FALSE}
*** This is the input file. All parameters that can be modified in the simulation
    are enclosed in this file. Variable names must be exact. To change the initial
    condition, edit the global.hpp file and run 'make'. ***

Select PDE ( "burgers" or "modified_burgers" or "advection" or "advection_diffusion" )
equation modified_burgers

Maximum number of wavelet levels in the simulation:
max_scale 8

Starting wavelet level (should be at least 2):
shift 2

Wavelet coefficient thresholding parameter
threshold 0.0000005

Half the number of interpolation points per stencil:
interp_points 2

Final time in the simulation:
tf 0.5

Number of timesteps in the simulation:
num_timesteps 8000

Advection velocity:
advec_vel 1.0

Coefficient of diffusivity:
diffusivity 0.01

Type of adjacent zone to add at each timestep
buffer_type Type1

Width of the buffer zone for Type1 wavelets
buffer_width 2

Height of the buffer zone for Type2 wavelets
buffer_height 1

Write data to file or not ( 1 = write, 0 = don't write ):
ifwrite 1
 ```

## Grid generation
The collocation points in the simulation are stored in a hierarchy of class objects
name `collPnt` and derived from the class `CollocationPoint`. The class is defined as:

```{r,eval=FALSE}
class CollocationPoint {

    //------- Public members -----------------------//
    public: 
        double x;                                   // location on the one-dimensional grid
        double u;                                   // the solution at the point
        double ux;                                  // the first derivative wrt x at the point
        double uxx;                                 // the second derivative wrt x at the point
        double scaling_coeff;                       // the scaling coefficient at the point
        double detail_coeff;                        // the wavelet coefficient at the point (if it is an odd point)
        bool isMask;                                // stores whether or not the point is in the mask at current iterate
        bool isOdd;                                 // stores whether or not the point corresponds to a wavelet
        bool isBuffer;                              // stores whether the point is in the adjacent zone
        bool isNew;                                 // determines whether corresponding wavelet needs to be computed
};
```
The grid points at each level $j = 0, \dots, J$ are computed as

\begin{equation}
x^{j}_{k} = 2^{-(j+\delta)} k
\end{equation}

for $k=0,\dots,2^{j+\delta}$, where $\delta$ is some integer shifting parameter, allowing one to dictate that the coarsest level of resolution have 
a smaller spacing than $1$.

## Wavelet transform & interpolation

## Wavelet construction

## Adaptive computational grid
Quantities of interest in the simulation are defined for each collocation point object. Among these are the state 
variable $u$, its first and second derivatives with respect to the spatial dimension $x$, $u_{x}$ and $u_{xx}$, as well
as the scaling coefficients $c_{k}^{j}$ and detail coefficients $d_{k}^{j}$. Each collocation point is assigned a boolean
value `isMask` which determines whether or not the point is a part of the adaptive computational grid at the current timestep.
Before the points are advanced in time, a buffer layer of grid points is added to the computational grid. These are adjacent 
wavelets which may become active during the next timestep. There are two strategies for adding the adjacent or buffer zone of 
wavelets: 

- `Type 1`, where points nearest an active wavelet are included. 
- `Type 2`, where insight about the physics of the problem is applied to predict where resolution will be required.
However the points must be labeled as `isBuffer`
## Timestepping

## Computation of spatial derivatives

## To-do list
- develop an implicit time scheme
- modify input file to pass left and right domain values
- modify input file to accept inhomogenous dirichlet data and periodic conditions
- Need to modify `write2file.cpp` to properly compute number of active points 
