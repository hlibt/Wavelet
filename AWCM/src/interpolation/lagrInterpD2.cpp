#include <math.h>
#include <cmath>
#include "../global.hpp"

double lagrInterpD2(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN) {
    double lagrDerivative=0.;
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
    for (int i=leftPnt;i<=rightPnt;i++) {
        double sumOuter=0.;
        for (int l=leftPnt;l<=rightPnt;l++) {
            if (l!=i) {
                double sumInner=0.;
                for (int m=leftPnt;m<=rightPnt;m++) {
                    if (m!=i && m!=l) {
                        double product=1.;
                        for (int k=leftPnt;k<=rightPnt;k++) {
                            if (k!=i && k!=l && k!=m) {
                                product*=(x-gridPnts[k])/(gridPnts[i]-gridPnts[k]);
                            }
                        }
                        sumInner+=product/(gridPnts[i]-gridPnts[m]);
                    }
                }
                sumOuter+=sumInner/(gridPnts[i]-gridPnts[l]);
            }
        }
        lagrDerivative+=funcPnts[i]*sumOuter;
    }
    return lagrDerivative;
}
