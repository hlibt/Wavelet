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
            gridPnts[j][k+N] = 2. * pow(2.,-(j+shift)) * k;             // x-locations of each collocation point
        }                                                       	    //
    }

    // generate scaling function   
    scaling_subd(phi,gridPnts,0,2,J,interpPnts);

    // compute the integrals of the scaling function at level j=0 using the trapezoid rule with lots of points
    double summation = 0.;
    for (int i=0;i<jPnts(J);i++) {
        if ( i == 0 || ( i == jPnts(J) - 1 ) ) {
            summation += phi[J][i] / 2.;
        } else {
            summation += phi[J][i];
        }
    }
    for (int i=0;i<jPnts(0);i++) {
        collPnt[0][i].integral = summation * abs( collPnt[J][1].x - collPnt[J][0].x );
    }

    // generate wavelets at all levels and then compute integral
    for (int j=0;j<J;j++) {
        detail_subd(psi,gridPnts,j,1,J,interpPnts);
        double summation = 0.;
        for (int i=0;i<jPnts(J);i++) {
            if ( i == 0 || ( i == jPnts(J) - 1 ) ) {
                summation += psi[J][i] / 2.;
            } else {
                summation += psi[J][i];
            }
        }
        for (int i=0;i<jPnts(j+1);i++) {
            if ( i%2==1 ) {
                collPnt[j+1][i].integral = summation * abs( collPnt[J][1].x - collPnt[J][0].x );
            }
        }
        for (int jstar=0;jstar<=J;jstar++) {
            for (int k=0;k<jPnts(jstar);k++) {
                psi[jstar][k] = 0.;
            }
        }
    }
        
}
