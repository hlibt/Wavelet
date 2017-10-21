#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../global.hpp"
using namespace std;

void thresholding(double** detCoeff,bool** mask,double epsilon,int Jmax) {
    for (int j=0;j<Jmax;j++) {                                          //
        int N=jPnts(j);                                                 //
        for (int k=0;k<N-1;k++) {                                       //
            if (abs(detCoeff[j][k])<epsilon) {                          //
                mask[j+1][2*k+1]=false;                                 // knock out points below threshold
                detCoeff[j][k]=0.;                                      // set to zeros ones not active
            } else {                                                    // 
                mask[j+1][2*k+1]=true;                                  // keep points above threshold, include in mask
            }                                                           //
        }                                                               //
    }                                                                   //
    return;                                                             //
}
