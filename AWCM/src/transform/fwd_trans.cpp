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
    double lagrange_coeff;                                              // the lagrange polynomial evaluated at interp point
    int n=pow(2,Jmax+1)+1;                                              //
    for (int i=0;i<n;i++) {                                             //
        c[Jmax][i]=u[i];                                                // set scaling coefficients to zero
    }                                                                   //
    for (int j=Jmax-1;j>=0;j--) {                                       // begin forward transform process
        int n=pow(2,j+1)+1;                                             // number of points at each level
        int Nstar;                                                      // adjustable copy of N
        if (j==0) {                                                     // 
            Nstar=1;                                                    // nearest points can only be 1
        } else if (j==1 && N>2) {                                       //
            Nstar=2;                                                    // cannot have nearest points greater than 2
        } else {                                                        //
            Nstar=N;                                                    // else allow number of points to be desired
        }                                                               //
        for (int k=0;k<n;k++) {                                       // 
            c[j][k]=c[j+1][2*k];                                        // even points stay the same
            int L1=-Nstar+1;                                            //
            int L2=Nstar;                                               //
            if (L1+k<0) {                                               //
                L2+=abs(0-(L1+k));                                      //
                L1+=abs(0-(L1+k));                                      //
            } else if (L2+k>=n) {                                       //
                L1-=abs((n-1)-(L2+k));                                  //
                L2-=abs((n-1)-(L2+k));                                  //
            }                                                           //
            double tmp=0.;                                              //
            for (int l=L1;l<=L2;l++) {                                  // 
                lagrange_coeff=lagrange_interp(x[j+1][2*k+1],       //
                                x[j],k+l,L1+k,L2+k);                //
             //   lagrange_coeff=0.5;
                tmp+=lagrange_coeff*c[j+1][2*k+2*l];                    //
            }                                                           //
            d[j][k]=.5*(c[j+1][2*k+1]-tmp);                             // detail coefficients
        }                                                               //
    }                                                                   //
    return;                                                             //
}                                                                       //
