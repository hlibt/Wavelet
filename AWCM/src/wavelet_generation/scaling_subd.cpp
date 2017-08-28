#include <iostream>
#include <cmath>
#include <math.h>
#include "wavelet_generation.hpp"
using namespace std;

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
        f[j][i]=kronecker_delta(i,m);                                       // set coefficients to at level j to kronecker delta function
    }                                                                       //
    for (int jstar=j;jstar<Jmax;jstar++) {                                  // begin inverse transform process started from level j
        int iPnts=setN(N,jstar);                                            // set number of interp points based on the level j
        int n=pow(2,jstar+1)+1;                                             // number of points at level jstar 
        for (int i=0;i<n-1;i++) {                                           // loop through all points but last in grid at level jstar
            double xEval=x[jstar+1][2*i+1];                                 // grid point for polynomial to be evaluated at
            f[jstar+1][2*i]=f[jstar][i];                                    // even points remain the same
            f[jstar+1][2*i+1]=lagrInterp(xEval,x[jstar],f[jstar],i,iPnts,n);// odd points
        }                                                                   //
        f[jstar+1][2*(n-1)]=f[jstar][n-1];                                  // last even point is the same
    }                                                                       //
    return;
}                                                                       

double lagrInterp(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN) {
    double sum=0.;
    int leftPnt=-n+1+i;
    int rightPnt=n+i;
    if ( leftPnt < 0 ) {
        rightPnt+=abs(0-leftPnt);
        leftPnt+=abs(0-leftPnt);
    } else if ( rightPnt > maxN-1 ) {
        leftPnt-=abs((maxN-1)-rightPnt);
        rightPnt-=abs((maxN-1)-rightPnt);
    }   
    for (int l=leftPnt;l<=rightPnt;l++) {
        double product=1.;
        for (int k=leftPnt;k<=rightPnt;k++) {
            if (k==l) {
            } else {
                product*=(x-gridPnts[k])/(gridPnts[l]-gridPnts[k]);
            }
        }
        sum+=product*funcPnts[l];
    }    
    return sum;
}

int setN(int n, int j) {
    if (j==0) {
        return 1;
    } else if (j==1 && n>2) {
        return 2;
    } else if (j==2 && n>3) {
        return 3;
    } else {
        return n;
    }
}

double kronecker_delta(int k, int m) {
    if (k==m) {
        return 1.;
    } else {
        return 0.;
    }
}   
