#include <iostream>
#include <math.h>
#include <cmath>

double lagrInterpD1(double x,double* gridPnts,double* funcPnts,int n) {
    double sum_outer=0.;
    for (int l=0;l<=n;l++) {
        double product=1.;
        double sum_inner=0.;
        for (int k=0;k<=n;k++) {
            if (k==l) {
            } else {
                product*=(x-gridPnts[k])/(gridPnts[l]-gridPnts[k]);
                sum_inner+=1./(x-gridPnts[k]);
            }
        }
        sum_outer+=product*sum_inner*funcPnts[l];
    }
    return sum_outer;
}

