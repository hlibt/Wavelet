#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include "../CollocationPoint.hpp"
#include "interpolation.hpp"
#include "../global.hpp"

vector<double> scaling_subd(vector<double> x,int j,int m) {
    
    //---------------------------------------------------------------------------//
    // Information: scaling_subd performs the interpolating subdivision algorithm
    //              in order to compute the scaling functions phi_j,m sampled at the 
    //              specific locations of x_Jmax, k. 
    //
    // Input: 
    //              j       - level of the scaling function
    //              m       - translaton parameter of the scaling function
    //              J       - maximum desired grid level for the point x (global variable)
    //              N       - half the number of nearest points to use in subdivision scheme (global variable)
    // Output:
    //              phi     - the scaling function of dilation j and translation k,
    //                        sampled at as many points as are grid level Jmax
    //--------------------------------------------------------------------------//     
                                                                                //
    vector<double> phi( jPnts(J) );                                             // create vector to store wavelet
    int n = jPnts(j);                                                           // number of points of phi at j level
    for (int i=0;i<n;i++) {                                                     // 
        phi[indexShift(J,j,i)] = kronecker_delta(i,m);                          // j coefficients to kronecker delta function
    }                                                                           //
    for (int jstar=j;jstar<Jmax;jstar++) {                                      // inverse transform process started from level j
        int n = jPnts(jstar);                                                   // number of points at level jstar 
        for (int i=0;i<n-1;i++) {                                               // loop through all points but last at J
            double xEval = x[jstar+1][2*i+1];                                   // grid point for polynomial to be evaluated at
            phi[jstar+1][2*i] = phi[jstar][i];                                  // even points remain the same
            phi[jstar+1][2*i+1] = lagrInterp(xEval,x[jstar],phi[jstar],i,N,n);  // odd points
        }                                                                       //
        phi[jstar+1][2*(n-1)] = phi[jstar][n-1];                                // last even point is the same
    }                                                                           //
    return;                                                                     //
}                                                                       

//------- Define the kronecker delta function ----------------------------------//
double kronecker_delta(int k, int m) {                                          //
    if (k==m) {                                                                 //
        return 1.;                                                              //
    } else {                                                                    //
        return 0.;                                                              //
    }                                                                           //
}                                                                               //
