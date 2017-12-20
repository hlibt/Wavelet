#include <iostream>
#include <cmath>
#include <math.h>
#include "../CollocationPoint.hpp"
#include "interpolation.hpp"
#include "../global.hpp"

void detail_subd(CollocationPoint** collPnt, double* detail_func,int j,int m) {
    
    //----------------------------------------------------------------------//
    // Information: detail_subd.cpp performs the interpolating subdivision algorithm
    //              in order to determine the wavelet psi_j,m sampled at the 
    //              specific locations of x_J,k. 
    //
    // Input: 
    //              collPnt - the collocation point objects
    //              scaling_func - the output scaling function
    //              j     - level of the scaling function
    //              m     - translaton parameter of the scaling function
    //              interpPnts - half the number of points in the interpolation stencil
    //
    // Output:
    //              phi_j,m(x_Jmax,k)
    //----------------------------------------------------------------------//     
 
    //------- Create temporary vector of collocation point objects ---------//
    CollocationPoint** cp = new CollocationPoint*[J+1];                     // copy the collocation point data 
    for (int jstar=0;jstar<=J;jstar++) {                                    //
        int n = jPnts(jstar);                                               //
        cp[jstar] = new CollocationPoint[n];                                //
    }                                                                       //

    //------- Copy x-locations and set scaling coefficients ----------------//
    for (int jstar=0;jstar<=J;jstar++) {                                    //
        int n = jPnts(jstar);                                               //
        for (int i=0;i<n;i++) {                                             //
            cp[jstar][i].x = collPnt[jstar][i].x;                           //
            cp[jstar][i].scaling_coeff = 0.;                                //
        }                                                                   //  
    }                                                                       //
    
    //------- Place unit impulse at level j index i=m ----------------------//
    for (int jstar=j;jstar<J;jstar++) {                                     // 
        int n = jPnts(jstar);                                               //
        for (int i=0;i<n-1;i++) {                                           //
            cp[jstar+1][2*i+1].detail_coeff = (jstar == j) * (i == m);      // set 'detail' coefficients
        }                                                                   //
    }                                                                       //

    //------- Build up wavelet to level J ----------------------------------//
    for (int jstar=j;jstar<J;jstar++) {                                     // begin inverse transform process
        int n = jPnts(jstar);                                               // number of points at level jstar
        for (int i=0;i<n-1;i++) {                                           // loop through all points but the last at level j
            double xeval = cp[jstar+1][2*i+1].x;                            // set the point to be evaluated at lagrange polynomial
            cp[jstar+1][2*i].scaling_coeff = cp[jstar][i].scaling_coeff;    // even points stay the same
            cp[jstar+1][2*i+1].scaling_coeff = 2. *                         //
                                    cp[jstar+1][2*i+1].detail_coeff +       // 
                                    lagrInterp(xeval,cp,jstar,i,n);         // 
        }                                                                   //
        cp[jstar+1][2*(n-1)].scaling_coeff = cp[jstar][n-1].scaling_coeff;  // last even point
    }                                                                       //

    //------- Transfer vector f[J] to scaling_func -------------------------//
    for (int i=0;i<jPnts(J);i++) detail_func[i] = cp[J][i].scaling_coeff;   // 

    //------- Cleanup ------------------------------------------------------//
    for (j=0;j<=J;j++) {
        delete[] cp[j];
    }
    delete[] cp;

    return;
}
