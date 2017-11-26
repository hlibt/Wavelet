#include <iostream>
#include <cmath>
#include <math.h>
#include "../global.hpp"
using namespace std;

double lagrInterp_old(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN) {
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
    if (leftPnt<0) {
        while ( leftPnt < 0 ) {
            n--;
            leftPnt=-n+1+i;
            rightPnt=n+i;
        }
        while ( rightPnt > (maxN-1) ) {
            n--;
            leftPnt=-n+1+i;
            rightPnt=n+1;
        } 
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
