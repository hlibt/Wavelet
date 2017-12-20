#include <math.h>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"

double lagrInterpD1(double x, CollocationPoint** collPnt, int j, int i, int numPnts) {
    
    double lagrDerivative = 0.;
    int leftPnt = -interpPnts + 1 + i;
    int rightPnt = interpPnts + i;
    while ( leftPnt < 0 ) {
        leftPnt++;
        rightPnt++;
    }
    while ( rightPnt > (numPnts-1) ) {
        leftPnt--;
        rightPnt--;
    }
    for (int p=leftPnt;p<=rightPnt;p++) {
        double sum = 0.;
        for (int k=leftPnt;k<=rightPnt;k++) {
            if (k!=p) {
                double productUpper = 1.;
                for (int l=leftPnt;l<=rightPnt;l++) {
                    if ( l!=k && l!=p ) { 
                        productUpper *= ( x - collPnt[j][l].x );
                    }
                }
                sum += productUpper;
            }
        }
        double productLower = 1.;
        for (int k=leftPnt;k<=rightPnt;k++) {
            if (k!=p) { 
                productLower *= ( collPnt[j][p].x - collPnt[j][k].x );
            }
        }
        lagrDerivative += collPnt[j][p].scaling_coeff * sum / productLower;
    }
    return lagrDerivative;
}
