#include <cmath>
#include "scaling_subd.hpp"

double* scaling_subd(int j,int m,int k,int Jmax,int npnts) {
    
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
    double weight=0.5;
    double** f=new double*[Jmax];                                       // function at points for interpolating subdivision
    for (int i=j;i<=Jmax;i++) {                                         //
        int n=pow(2,i+1)+1;                                             // number of points at level j
        f[i]=new double[n];                                             // intialize columns of f
    }                                                                   //
    int n=pow(2,j+1)+1;                                                 //
    for (int i=0;i<n;i++) {                                             //
        f[j][i]=kronecker_delta(k,m);                                   // set coefficients to kronicker delta function
    }                                                                   //
    for (int jstar=j;jstar<=Jmax;jstar++) {                             // begin inverse transform process
        int n=pow(2,jstar+1)+1;                                         // number of points at level jstar
        for (int i=0;i<n;i++) {                                         // 
            f[jstar+1][2*i]=f[jstar][i];                                // even points stay the same
            double tmp=0.;                                              // summation variable
            for (int l=-npnts+1;l<=npnts;l++) {                         // 
                tmp+=weight*f[jstar][i+l];                              //
            }                                                           //
            f[jstar+1][2*i+1]=tmp;                                      // odd points
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
