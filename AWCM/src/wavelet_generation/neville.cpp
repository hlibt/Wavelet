#include <iostream>
#include <cmath>
#include <math.h>
#include "wavelet_generation.hpp"
using namespace std;

double neville(double y,double* x,double* f,int n) {
    double** Q=new double*[n];
    for (int i=0;i<n;i++) Q[i]=new double[n];
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            if (j==0) {
                Q[i][j]=f[i];
            } else {
                Q[i][j]=0.;
            }
        }
    }
    for (int i=1;i<n;i++) {
        for (int j=1;j<=i;j++) {
            Q[i][j]=( (y-x[i-j])*Q[i][j-1] - (y-x[i])*Q[i-1][j-1] ) / (x[i]-x[i-j]);
        }
    }
    return Q[n-1][n-1];
}
    
