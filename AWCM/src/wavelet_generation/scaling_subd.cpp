#include <iostream>
#include <cmath>
#include <math.h>
#include "wavelet_generation.hpp"
using namespace std;

double* scaling_subd(double** x,int j,int m,int Jmax,int N) {
    
    //------------------------------------------------------------------//
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
    //------------------------------------------------------------------//     
                                                                        //
    double w;                                                           // the lagrange polynomial evaluated at interp point
    double** f=new double*[Jmax];                                       // function at points for interpolating subdivision
    for (int i=0;i<=Jmax;i++) {                                         //
        int n=pow(2,i+1)+1;                                             // number of points at level j
        f[i]=new double[n];                                             // intialize columns of f
    }                                                                   //
    int n=pow(2,j+1)+1;                                                 //
    for (int i=0;i<n;i++) {                                             //
        f[j][i]=kronecker_delta(i,m);                                   // set coefficients to kronecker delta function
    }                                                                   //
    for (int jstar=j;jstar<Jmax;jstar++) {                              // begin inverse transform process
        int n=pow(2,jstar+1)+1;                                         // number of points at level jstar 
        for (int i=0;i<n;i++) {                                         // loop through all points in grid at level jstar
            f[jstar+1][2*i]=f[jstar][i];                                // even points remain the same
            double sum=0.;                                              // summation variable
            for (int l=-N+1;l<=N;l++) {                                 // 
                if ( l+i >=0 && l+i < n ) {                             //
                    double c=x[jstar+1][2*i+1];                         // constant
                    w=lagrange_interp(c,x[jstar],i+l,-N+1+i,N+i,n);     //
                    cout << w << endl;
                    sum+=w*f[jstar][l+i];                               //
                }
            }                                                           //
            f[jstar+1][2*i+1]=sum;                                      // odd points
        }                                                               //
    }                                                                   //
    return f[Jmax];                                                     // the final scaling function at sampled points
}                                                                       //

double lagrange_interp(double eval_point,double* x,int i,int N1,int N2,int n) {
    double prod=1.;
    for (int k=N1;k<=N2;k++) {
        if ( k>=0 && k<n ) {
            if (k==i) {
            } else {
                prod*=(eval_point-x[k])/(x[i]-x[k]);
            }
        }
    }
    return prod;
}

double kronecker_delta(int k, int m) {
    if (k==m) {
        return 1.;
    } else {
        return 0.;
    }
}
