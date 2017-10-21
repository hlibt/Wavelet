#include <iostream>
#include <cmath>
#include <math.h>
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
using namespace std;

void lftd_fwd_trans(double** x,double** scalCoeff,double** detCoeff,int Jmax,int N_even,int N_odd) {
    
    //------------------------------------------------------------------//
    // Information: fwd_trans performs the forward wavelet 
    //              transformation to recover the scaling and detail
    //              coefficients at each lower level until j=0. 
    //
    // Input: 
    //              x     - grid at each level j
    //              u     - solution at finest resolution Jmax
    //              scalCoeff - scaling coeff's 
    //              detCoeff - detail coeff's
    //              Jmax  - maximum grid level
    //              N     - half the number of nearest neighbors to use
    //------------------------------------------------------------------//     
    
    for (int j=Jmax-1;j>=0;j--) {                                       // begin forward transform process
        int n=jPnts(j);                                                 // number of points at level j
        double xEval;                                                   // interpolant evaluation point
        double* copy = new double[n];                                   // copy even scaling coeff's from level j+1
        for (int i=0;i<n;i++) copy[i]=scalCoeff[j+1][2*i];              // copy
        for (int i=0;i<n-1;i++) {                                       //
            xEval=x[j+1][2*i+1];                                        // point for polynomials to be evaluated 
            detCoeff[j][i]=.5*(scalCoeff[j+1][2*i+1]-lagrInterp(xEval,  // detail coefficient computation
                        x[j],copy,i,N_even,n));                         //
        }                                                               //
        for (int i=0;i<n-1;i++) {                                       // 
            xEval=x[j][k];                                              // even point for evaluation
            scalCoeff[j][i]=scalCoeff[j+1][2*i]+lagrInterp(xEval,x[j],  // determine scaling coefficients
                        detCoeff[j],i,N_odd,n-1);   
        }                                                               //
        xEval=x[j][n-1];                                                // evaluation point for last scaling coefficient
        scalCoeff[j][n-1]=scalCoeff[j+1][2*(n-1)]+lagrInterp(xEval,     // compute last scaling coefficient
                        x[j],detCoeff[j],n-2,N_odd,n-1);                //
        delete[] copy;                                                  // cleanup memory
    }                                                                   //
    return;                                                             //
}  
