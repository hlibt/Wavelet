#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../global.hpp"
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
using namespace std;

void wavelet_integrals(CollocationPoint** collPnt) {

    // create arrays
    double** gridPnts = new double*[J+1];                          
    double** phi = new double*[J+1];
    double** psi = new double*[J+1];
    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        gridPnts[j] = new double[N];
        phi[j] = new double[N];
        psi[j] = new double[N];
    }
    for (int j=0;j<=J;j++) {                                   
        int N = ( jPnts(j) - 1 ) / 2;
        for (int k=-N;k<=N;k++) {                                      	//
            gridPnts[j][k+N] = 0.5 + pow(2.,-(j+shift)) * k;             // x-locations of each collocation point
        }                                                       	    //
    }

    // compute the integrals of the scaling function at level j=0 using the trapezoid rule with lots of points
    for (int l=0;l<jPnts(0);l++) {
        scaling_subd(phi,gridPnts,0,l,J,interpPnts);                    // generate the scaling function
        double summation = 0.;                                          // summation variable for trapezoid rule
        for (int i=0;i<jPnts(J);i++) {
            if ( i == 0 || ( i == jPnts(J) - 1 ) ) {
                summation += phi[J][i] / 2.;
            } else {
                summation += phi[J][i];
            }
        }
        collPnt[0][l].integral = summation * abs( gridPnts[J][1] - gridPnts[J][0] );
    }
    
    // generate wavelets at all levels and then compute integral
    for (int j=0;j<J;j++) {
        int N = jPnts(j);
        for (int l=0;l<N;l++) {
            if ( l < N -1 && collPnt[j+1][2*l+1].isMask == true ) {
                detail_subd(psi,gridPnts,j,l,J,interpPnts);
                double summation = 0.;                                      // summation variable for trapezoid rule
                for (int i=0;i<jPnts(J);i++) {
                    if ( i == 0 || ( i == jPnts(J) - 1 ) ) {
                        summation += psi[J][i] / 2.;
                    } else {
                        summation += psi[J][i];
                    }
                }
                collPnt[j+1][2*l+1].integral = summation * abs( gridPnts[J][1] - gridPnts[J][0] );
            }
        }
    }
        
}

