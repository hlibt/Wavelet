#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../global.hpp"
#include "../CollocationPoint.hpp"
using namespace std;

double integrate(CollocationPoint** collPnt) {

    // summation variable to compute the final integral of the function
    double summation = 0.;

    // sum up the integrals of the scaling functions (j=0) multiplied by their corresponding weights
    for (int i=0;i<jPnts(0);i++) {
        summation += collPnt[0][i].integral * collPnt[0][i].scaling_coeff;
    }

        // sum up integrals from the detail wavelets multiplied by their detail coefficients (at all levels)
    for (int j=1;j<=J;j++) {
        for (int i=0;i<jPnts(j);i++) {
            if ( collPnt[j][i].isMask == true ) {
               summation += collPnt[j][i].integral * collPnt[j][i].detail_coeff;
            }
        }
    }

    return summation; 
}
