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

    // sum up the integrals of the scaling wavelets at level j=0
    for (int i=0;i<jPnts(0);i++) {
        summation += collPnt[0][i].integral;
    }

    // sum up integrals from the detail wavelets at all levels
    for (int j=1;j<=J;j++) {
        int n = jPnts(j);
        for (int i=0;i<n;i++) {
            if ( collPnt[j][i].isMask == true ) {
               summation += collPnt[j][i].integral;
            }
        }
    }
        
    return summation; 
}
