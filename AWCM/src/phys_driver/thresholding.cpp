#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void thresholding(CollocationPoint** collPnt, double epsilon) {
    for (int i=0;i<jPnts(0);i++) {                                          //
        collPnt[0][i].isMask = true;                                        //
    }                                                                       //
    for (int j=1;j<=J;j++) {                                                //
        int N = jPnts(j);                                                   //
        for (int k=0;k<N;k++) {                                             //
            if ( k%2==1 && abs(collPnt[j][k].detail_coeff) > epsilon ) {    // if point is odd and if the detail coefficient is large enough
                collPnt[j][k].isMask = true;                                //
            }                                                               //
        }                                                                   //
    }                                                                       //
    return;                                                                 //
}
