#include <iostream>
#include <cmath>
#include <math.h>
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
using namespace std;

void fwd_trans(double** x,double** scalCoeff,double** detCoeff,int Jmax,int N) {
    
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
        for (int i=0;i<n;i++) {                                         // 
            scalCoeff[j][i]=scalCoeff[j+1][2*i];                        // even scaling coeff's the same
        }                                                               //
        for (int i=0;i<n-1;i++) {                                       //
            double xEval=x[j+1][2*i+1];                                 // point for polynomials to be evaluated 
            detCoeff[j][i]=.5*(scalCoeff[j+1][2*i+1]-                   // set detail coefficients
                        lagrInterp(xEval,x[j],scalCoeff[j],i,N,n));     //
        }                                                               //
    }                                                                   //
    return;                                                             //
}  
