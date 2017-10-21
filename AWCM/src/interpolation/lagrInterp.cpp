#include <iostream>
#include <cmath>
#include <math.h>
#include "../global.hpp"
using namespace std;

double lagrInterp(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN) {
    double sum=0.;
    int leftPnt=-n+1+i;
    int rightPnt=n+i;
    while ( leftPnt < 0 ) {
        leftPnt++;
        rightPnt++;
    }
    while ( rightPnt > (maxN-1) ) {
        leftPnt--;
        rightPnt--;
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
