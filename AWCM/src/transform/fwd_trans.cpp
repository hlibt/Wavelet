#include <iostream>
#include <cmath>
#include <math.h>
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
#include "../CollocationPoint.hpp"
using namespace std;

void fwd_trans(CollocationPoint** collPnt,double* inputSig,int Jmax,int N) {
    
    //--------------------------------------------------------------------------//
    // Information: fwd_trans performs the forward wavelet                  
    //              transformation to compute the scaling and detail
    //              coefficients for all levels 
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              inputSig                - the input function to be decomposed
    //              Jmax                    - maximum grid level
    //              N                       - half the number of interp. points
    //--------------------------------------------------------------------------//     
    
    int n = jPnts(Jmax);                                                        // number of points at level Jmax
    for (int i=0;i<n;i++) {                                                     // loop through the points at level Jmax
        collPnt[Jmax][i].scaling_coeff = inputSig[i];                           // set the scaling coefficients at Jmax to the input
    };                                                                          //
    for (int j=Jmax-1;j>=0;j--) {                                               // copy even coefficients from level j+1 to j
        int n = jPnts(j);                                                       // 
        for (int i=0;i<n;i++) {                                                 // loop through points at level j
            collPnt[j][i].scaling_coeff = collPnt[j+1][2*i].scaling_coeff;      // even scaling coeff's the same, just copy them
        }                                                                       //
        for (int i=0;i<n-1;i++) {                                               //
            double xeval = collPnt[j+1][2*i+1].x;                               // define the point for polynomials to be evaluated 
            collPnt[j][i].detail_coeff = .5*( collPnt[j+1][2*i+1].scaling_coeff // compute detail coefficients
                        - lagrInterp(xeval,collPnt[j],i,N,n) );                 //
        }                                                                       //
    }                                                                           //
    return;                                                                     //
}  
