#include <iostream>
#include <cmath>
#include <math.h>
#include "scaling_subd.hpp"
using namespace std;

double* scaling_subd(double** x,int j,int m,int Jmax,int npnts) {
    
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
    //              npnts - half the number of nearest points to use in subdivision scheme
    // Output:
    //              phi_j,m(x_Jmax,k)
    //------------------------------------------------------------------//     
                                                                        //
    double lagrange_coeff;                                              // the lagrange polynomial evaluated at interp point
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
        for (int i=0;i<n-1;i++) {                                       // 
            f[jstar+1][2*i]=f[jstar][i];                                // even points stay the same
            double tmp=0.;                                              // summation variable
            int L1=-npnts+1;
            int L2=npnts;
            if (L1+i<0) {                                               //
                L1+=abs(0-L1+i);                                        //
                L2+=abs(0-L1+i);                                        //
            }                                                           //
            if (L2+i>=n) {                                              //
                L2-=abs((n-1)-(L2+i));                                  //
                L1-=abs((n-1)-(L2+i));                                  //
            }                                                           //
            for (int l=L1;l<=L2;l++) {                                  // 
                lagrange_coeff=lagrange_interp(x[jstar+1][2*i+1],       //
                                x[jstar],l+i,L1,L2);                    //
                tmp+=lagrange_coeff*f[jstar][i+l];                      //
            }                                                           //
            f[jstar+1][2*i+1]=tmp;                                      // odd points
        }                                                               //
    }                                                                   //
    return f[Jmax];                                                     // the final scaling function at sampled points
}                                                                       //

double lagrange_interp(double eval_point,double* x,int i,int N1,int N2) {
    double prod=1.;
    for (int k=N1;k<=N2;k++) {
        if (k==i) {
        } else {
            prod*=(eval_point-x[k])/(x[i]-x[k]);
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
