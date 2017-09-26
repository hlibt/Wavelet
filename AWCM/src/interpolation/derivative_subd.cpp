#include <iostream>
#include <cmath>
#include <math.h>
#include "wavelet_generation.hpp"
#include "../interpolation/interpolation.hpp"

void derivative_subd(double** df,double** f,double** x,int j,int N) {
    
    //----------------------------------------------------------------------//
    // Information: scaling_subd performs the interpolating subdivision algorithm
    //              in order to determine the scaling functions phi_j,m sampled at the 
    //              specific locations of x_Jmax,k. 
    //
    // Input: 
    //              f - the scaling coefficients   
    //              j     - level of the grid points
    //              Jmax  - maximum desired grid level for the point x
    //              N - half the number of nearest points to use in subdivision scheme
    // Output:
    //----------------------------------------------------------------------//     
                                                                            //
    int n=pow(2,jstar+1)+1;                                                 // number of points at level jstar 
    for (int i=0;i<n-1;i++) {                                               // loop through all points but level jstar
        double xEval=x[jstar+1][2*i+1];                                     // grid point for polynomial to be evaluated at
        f[jstar+1][2*i+1]=lagrInterpD1(xEval,x[jstar],f[jstar],i,N,n);      // odd points
    }                                                                       //
    return;
}                                                                       
