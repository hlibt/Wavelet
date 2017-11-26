#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../global.hpp"
#include "../CollocationPoint.hpp"
using namespace std;

void thresholding(CollocationPoint** collPnt,double epsilon) {
    for (int j=J;j>0;j--) {                                             // loop through each level, starting from finest
        int n = jPnts(j);                                               // number of points at level j
        for (int k=0;k<n;k++) {                                         // loop through detail coefficients in each row
            if ( k%2 == 1 && collPnt[j][k].detail_coeff > epsilon ) {   // for the odd points, determine if the detail coeffcient is significant
                collPnt[j][k].isMask = true;                            // if so, set 'isMask' to true
            }                                                           //
        }                                                               //
    }                                                                   //
    return;                                                             // 
}
