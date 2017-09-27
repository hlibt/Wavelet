#include <iostream>
#include <cmath>
#include <math.h>
#include "interpolation.hpp"

void scaling_subd(double** f,double** x,int j,int m,int Jmax,int N) {
    
    //----------------------------------------------------------------------//
    // Information: scaling_subd performs the interpolating subdivision algorithm
    //              in order to determine the scaling functions phi_j,m sampled at the 
    //              specific locations of x_Jmax,k. 
    //
    // Input: 
    //              j     - level of the scaling function
    //              m     - translaton parameter of the scaling function
    //              Jmax  - maximum desired grid level for the point x
    //              k     - spatial index of x
    //              N - half the number of nearest points to use in subdivision scheme
    // Output:
    //              phi_j,m(x_Jmax,k)
    //----------------------------------------------------------------------//     
                                                                            //
    int n=pow(2,j+1)+1;                                                     // number of points of phi at j level
    for (int i=0;i<n;i++) {                                                 // 
        f[j][i]=kronecker_delta(i,m);                                       // j coefficients to kronecker delta function
    }                                                                       //
    for (int jstar=j;jstar<Jmax;jstar++) {                                  // inverse transform process started from level j
        int n=pow(2,jstar+1)+1;                                             // number of points at level jstar 
        for (int i=0;i<n-1;i++) {                                           // loop through all points but last at J
            double xEval=x[jstar+1][2*i+1];                                 // grid point for polynomial to be evaluated at
            f[jstar+1][2*i]=f[jstar][i];                                    // even points remain the same
            f[jstar+1][2*i+1]=lagrInterp(xEval,x[jstar],f[jstar],i,N,n);    // odd points
        }                                                                   //
        f[jstar+1][2*(n-1)]=f[jstar][n-1];                                  // last even point is the same
    }                                                                       //
    return;
}                                                                       

double kronecker_delta(int k, int m) {
    if (k==m) {
        return 1.;
    } else {
        return 0.;
    }
}   
