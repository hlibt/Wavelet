#include <iostream>
#include <cmath>
#include <math.h>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
#include "interpolation.hpp"

void scaling_subd(CollocationPoint** collPnt, double* scaling_func,int j,int m) {
    
    //----------------------------------------------------------------------//
    // Information: scaling_subd performs the interpolating subdivision algorithm
    //              in order to determine the scaling functions phi_j,m sampled at the 
    //              specific locations of x_Jmax,k. 
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

    //------- Copy x-locations ---------------------------------------------//
    for (int jstar=0;jstar<=J;jstar++) {                                    //
        int n = jPnts(jstar);                                               //
        for (int i=0;i<n;i++) {                                             //
            cp[jstar][i].x = collPnt[jstar][i].x;                           //
        }                                                                   //  
    }                                                                       //

    //------- Place unit impulse at level j index i=m ----------------------//
    int n = jPnts(j);                                                       // number of points of phi at j level
    for (int i=0;i<n;i++) {                                                 // 
        cp[j][i].scaling_coeff = ( i == m );                                // j coefficients to kronecker delta function
    }                                                                       //

    //------- Build up scaling function to level J -------------------------//
    for (int jstar=j;jstar<J;jstar++) {                                     // inverse transform process started from level j
        int n = jPnts(jstar);                                               // number of points at level jstar 
        for (int i=0;i<n-1;i++) {                                           // loop through all points but last at J
            double xeval = cp[jstar+1][2*i+1].x;                            // grid point for polynomial to be evaluated at
            cp[jstar+1][2*i].scaling_coeff = cp[jstar][i].scaling_coeff;    // even points remain the same
            cp[jstar+1][2*i+1].scaling_coeff = lagrInterp(xeval,cp,jstar,   // odd points
                                                    i,n);                   // 
        }                                                                   //
        cp[jstar+1][2*(n-1)].scaling_coeff = cp[jstar][n-1].scaling_coeff;  // last even point is the same
    }                                                                       //

    //------- Transfer vector f[J] to scaling_func -------------------------//
    for (int i=0;i<jPnts(J);i++) scaling_func[i] = cp[J][i].scaling_coeff;  // 

    //------- Cleanup ------------------------------------------------------//
    for (j=0;j<=J;j++) {
        delete[] cp[j];
    }
    delete[] cp;

    return;
}
