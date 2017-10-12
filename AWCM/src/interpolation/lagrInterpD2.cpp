#include <math.h>
#include <cmath>
#include "../global.hpp"

double lagrInterpD1(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN) {
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
/*    for (int l=leftPnt;l<=rightPnt;l++) {
        double product=1.;
        for (int k=leftPnt;k<=rightPnt;k++) {if (k!=l) product*=(x-gridPnts[k])/(gridPnts[l]-gridPnts[k]);}
        double sum=0.;
        for (int k=leftPnt;k<=rightPnt;k++) {if (k!=l) sum+=1./(x-gridPnts[k]);}
        lagrDerivative+=product*sum*funcPnts[l];
    } */
    for (int j=leftPnt;j<=rightPnt;j++) {
        double sum=0.;
        for (int k=leftPnt;k<=rightPnt;k++) {
            if (k!=j) {
                double productUpper=1.;
                for (int l=leftPnt;l<=rightPnt;l++) {if (l!=k && l!=j) productUpper*=(x-gridPnts[l]);}
                sum+=productUpper;
            }
        }
        double productLower=1.;
        for (int k=leftPnt;k<=rightPnt;k++) {if (k!=j) productLower*=(gridPnts[j]-gridPnts[k]);}
        lagrDerivative+=funcPnts[j]*sum/productLower;
    }
    return lagrDerivative;
}
