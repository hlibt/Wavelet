#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../global.hpp"
#include "../CollocationPoint.hpp"
using namespace std;

void thresholding(CollocationPoint** collPnt,double epsilon) {
    for (int j=J;j>0;j--) {                                             //
        int n = jPnts(j);                                               //
        for (int k=0;k<n;k++) {                                         // loop through detail coefficients in each row
            if ( k%2 == 1 && collPnt[j][k].detail_coeff < epsilon ) {   //
                collPnt[j][k].isMask = true;                            //
            }                                                           //
        }                                                               //
    }                                                                   //
    return;                                                             //
}
