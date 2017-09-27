#include <iostream>
#include <cmath>
#include <math.h>
#include "interpolation.hpp"

void detail_subd(double** f,double** x,int j,int m,int Jmax,int N) {
    
    //------------------------------------------------------------------//
    // Information: detail_subd performs the interpolating subdivision algorithm
    //              in order to determine the daughter functions psi_j,m sampled at the 
    //              specific locations of x_Jmax,k. 
    //
    // Input: 
    //              j     - level of the scaling function
    //              m     - translaton parameter of the scaling function
    //              Jmax  - maximum desired grid level for the point x
    //              N     - half the number of nearest points to use in subdivision scheme
    // Output:
    //              psi_j,m(x_Jmax,k)
    //------------------------------------------------------------------//     
                                                                        //
    double** d=new double*[Jmax+1];                                     // detail function coefficients
    for (int i=0;i<=Jmax;i++) {                                         //
        int n=pow(2,i+1)+1;                                             // number of points at level j
        d[i]=new double[n];                                             // initialize columns of d
    }                                                                   //
    int n=pow(2,j+1)+1;                                                 //
    for (int i=0;i<n;i++) {                                             //
        f[j][i]=0;                                                      // set 'scaling' coefficients to zero
    }                                                                   //
    for (int jstar=j;jstar<Jmax;jstar++) {                              // 
        int n=pow(2,jstar+1)+1;                                         //
        for (int i=0;i<n;i++) {                                         //
            d[jstar][i]=kronecker_delta(jstar,j)*kronecker_delta(i,m);  // set 'detail' coefficients
        }                                                               //
    }                                                                   //
    for (int jstar=j;jstar<Jmax;jstar++) {                              // begin inverse transform process
        int n=pow(2,jstar+1)+1;                                         // number of points at level jstar
        for (int i=0;i<n-1;i++) {                                       // loop through all points but the last at level j
            double xEval=x[jstar+1][2*i+1];                             // set the point to be evaluated at lagrange polynomial
            f[jstar+1][2*i]=f[jstar][i];                                // even points stay the same
            f[jstar+1][2*i+1]=2.*d[jstar][i]+lagrInterp(xEval,x[jstar], // odd points
                                f[jstar],i,N,n);                        // 
        }                                                               //
        f[jstar+1][2*(n-1)]=f[jstar][n-1];                              // last even point
    }                                                                   //
    delete[] d;                                                         // delete dynamic memory
    return;                                                             // 
}                                                                       
