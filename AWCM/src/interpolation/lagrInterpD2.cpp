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
    for (int j=leftPnt;j<=rightPnt;j++) {
        double sum=0.;
        for (int k=leftPnt;k<=rightPnt;k++) {
            if (k!=j) {
                double sum1=0.;
                double sum2=0.;
                for (int l=leftPnt;l<=rightPnt;l++) {
                    if (l!=k && l!=j) {
                        double productUpper=1.;
                        for (int lstar=leftPnt;lstar<=rightPnt;lstar++) {
                            if (lstar!=j && lstar!=l) {
                                productUpper*=gridPnts[l]*gridPnts[lstar];
                            }    
                        }
                        sum1+=gridPnts[l];
                        sum2+=productUpper;
                    }
                }
                sum+=3.*pow(x,2.)-2.*x*sum1+sum2;
            }
        }
        double productLower=1.;
        for (int k=leftPnt;k<=rightPnt;k++) {if (k!=j) productLower*=(gridPnts[j]-gridPnts[k]);}
        lagrDerivative+=funcPnts[j]*sum/productLower;
    }
    return lagrDerivative;
}
