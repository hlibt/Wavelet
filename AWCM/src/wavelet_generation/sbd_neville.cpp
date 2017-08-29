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
        for (int i=0;i<n;i++) {                                         //
            f[jstar+1][2*i]=f[jstar][i];                                // even points remain the same
            f[jstar+1][2*i+1]=neville(x[jstar+1][2*i+1],x[jstar],       // interpolate to obtain odd points
                                f[jstar],n);                            //
        }                                                               //
    }                                                                   //
    return f[Jmax];                                                     // the final scaling function at sampled points
}                                                                       //

double kronecker_delta(int k, int m) {
    if (k==m) {
        return 1.;
    } else {
        return 0.;
    }
}

