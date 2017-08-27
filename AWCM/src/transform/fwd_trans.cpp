#include <iostream>
#include <cmath>
#include <math.h>
#include "../wavelet_generation/wavelet_generation.hpp"
using namespace std;

void fwd_trans(double** x,double* u,double** c,double** d,int Jmax,int N) {
    
    //------------------------------------------------------------------//
    // Information: fwd_trans performs the forward wavelet 
    //              transformation to recover the scaling and detail
    //              coefficients at each lower level until j=0. 
    //
    // Input: 
    //              x     - grid at each level j
    //              u     - solution at finest resolution Jmax
    //              c     - scaling coeff's 
    //              d     - detail coeff's
    //              Jmax  - maximum grid level
    //              N     - half the number of nearest neighbors to use
    //------------------------------------------------------------------//     
                                                                        //
    int n=pow(2,Jmax+1)+1;                                              //
    for (int i=0;i<n;i++) {                                             //
        c[Jmax][i]=u[i];                                                // set scaling coefficients to zero
    }                                                                   //
    for (int j=Jmax-1;j>=0;j--) {                                       // begin forward transform process
        int n=pow(2,j+1)+1;                                             // number of points at level j
        int iPnts=setN(N,j);                                            // set number of interp points based on level j
        for (int i=0;i<n;i++) {                                         // 
            c[j][i]=c[j+1][2*i];                                        // even scaling coeff's the same
        }                                                               //
        for (int i=0;i<n-1;i++) {                                       //
            double xEval=x[j+1][2*i+1];                                 // point for polynomials to be evaluated 
            d[j][i]=.5*(c[j+1][2*i+1]-lagrInterp(xEval,x[j],c[j],       // calculate detail coefficients
                     i,iPnts,n));                                       //
        }                                                               //
    }                                                                   //
    return;                                                             //
}  
