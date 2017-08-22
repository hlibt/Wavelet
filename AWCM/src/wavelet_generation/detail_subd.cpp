#include <iostream>
#include <cmath>
#include <math.h>
#include "wavelet_generation.hpp"
using namespace std;

double* detail_subd(double** x,int j,int m,int Jmax,int N) {
    
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
    double** c=new double*[Jmax];                                       // scaling function coefficients
    double** d=new double*[Jmax];                                       // detail function coefficients
    for (int i=0;i<=Jmax;i++) {                                         //
        int n=pow(2,i+1)+1;                                             // number of points at level j
        c[i]=new double[n];                                             // intialize columns of c
        d[i]=new double[n];                                             // initialize columns of d
    }                                                                   //
    int n=pow(2,j+1)+1;                                                 //
    for (int i=0;i<n;i++) {                                             //
        c[j][i]=0;                                                      // set scaling coefficients to zero
    }                                                                   //
    for (int jstar=j;jstar<Jmax;jstar++) {                              //
        int n=pow(2,jstar+1)+1;                                         //
        for (int i=0;i<n;i++) {                                         //
            d[jstar][i]=kronecker_delta(jstar,j)*kronecker_delta(i,m);  //
        }                                                               //
    }                                                                   //
    for (int jstar=j;jstar<Jmax;jstar++) {                              // begin inverse transform process
        int n=pow(2,jstar+1)+1;                                         // number of points at level jstar
        int Nstar;                                                      // adjustable copy of N
        if (jstar==0) {                                                 // 
            Nstar=1;                                                    // nearest points can only be 1
        } else if (jstar==1 && N>2) {                                   //
            Nstar=2;                                                    // cannot have npnts greater than 2
        } else {                                                        //
            Nstar=N;                                                    //
        }                                                               //
        for (int i=0;i<n;i++) {                                         // 
            c[jstar+1][2*i]=c[jstar][i];                                // even points stay the same
            int L1=-Nstar+1;                                            //
            int L2=Nstar;                                               //
            if (L1+i<0) {                                               //
                L2+=abs(0-(L1+i));                                      //
                L1+=abs(0-(L1+i));                                      //
            } else if (L2+i>=n) {                                       //
                L1-=abs((n-1)-(L2+i));                                  //
                L2-=abs((n-1)-(L2+i));                                  //
            }                                                           //
            double tmp=0.;                                              //
            for (int l=L1;l<=L2;l++) {                                  // 
                lagrange_coeff=lagrange_interp(x[jstar+1][2*i+1],       //
                                x[jstar],i+l,L1+i,L2+i);                //
                tmp+=lagrange_coeff*c[jstar][i+l];                      //
            }                                                           //
            c[jstar+1][2*i+1]=2*d[jstar][i]+tmp;                        // odd points
        }                                                               //
    }                                                                   //
    return c[Jmax];                                                     // the final scaling function at sampled points
}                                                                       //
