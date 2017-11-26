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

    // number of trapezoid rule integration 
    int Jtrap = 10;
    int Ntrap = jPnts( Jtrap );

    // create arrays
    double** gridPnts = new double*[Jtrap+1];                          
    double** phi = new double*[Jtrap+1];
    double** psi = new double*[Jtrap+1];
    for (int j=0;j<=Jtrap;j++) {
        int N = jPnts(j);
        gridPnts[j] = new double[N];
        phi[j] = new double[N];
        psi[j] = new double[N];
    }
    for (int j=0;j<=J;j++) {                                   
        int N = jPnts(j);                                               
        for (int k=-N;k<=N;k++) {                                      	//
            gridPnts[j][k+N] = 2. * pow(2.,-(j+shift)) * k;             // x-locations of each collocation point
        }                                                       	    //
    }

    // generate scaling function   
    scaling_subd(phi,gridPnts,0,2,Jtrap,interpPnts);

    // compute the integrals of the scaling function at level j=0 using the trapezoid rule with lots of points
    double summation = 0.;
    for (int i=0;i<jPnts( Jtrap );i++) {
        if ( i == 0 || ( i == jPnts( Jtrap ) - 1 ) ) {
            summation += phi[Jtrap][i] / 2;
        } else {
            summation += phi[Jtrap][i];
        }
        cout << phi[Jtrap][i] << endl;
    }
    for (int i=0;i<jPnts(0);i++) {
        collPnt[0][i].integral = summation;
    }
    cout << summation << endl;

    // generate wavelets at all levels and then compute integral
    for (int j=0;j<J;j++) {
        detail_subd(psi,gridPnts,j,3,Jtrap,interpPnts);
        double summation = 0.;
        for (int i=0;i<jPnts( Jtrap );i++) {
            if ( i == 0 || ( i == jPnts( Jtrap ) - 1 ) ) {
                summation += psi[Jtrap][i] / 2.;
            } else {
                summation += psi[Jtrap][i];
            }
        }
        for (int i=0;i<jPnts(j+1);i++) {
            if ( i%2==1 ) {
                collPnt[j+1][i].integral = summation;
            }
        }
    }
        
}