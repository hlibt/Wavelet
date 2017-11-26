#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "global.hpp"
#include "CollocationPoint.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
using namespace std;

    //------------------------------------------------------------------------------//
    //                                                                              //
    //     This program tests wavelet collocation as an adaptive quadrature         //
    //     method for integrating nearly singular integrands                        //
    //                                                                              //
    //     AUTHOR: BRANDON GUSTO                                                    //
    //     DATE CREATED: Nov 21, 2017                                               //
    //                                                                              //
    //------------------------------------------------------------------------------//

int shift;
int J;
int interpPnts;

int main(void) {

    //------- Grid and tolerance parameters ----------------------------//
    shift = 2;                                                          // increases number of points on level j=0 (global variable)
    J = 7;                                                              // number of scales in the system
    interpPnts = 2;                                                     // half the number of points used for interpolation (2*interpPnts + 1)
    double threshold = pow(10.,-6);                                	    // error tolerance for wavelet coefficients (determines accuracy of solution)
    int i;                                                              // the usual counter variable for spatial index
    int j;                                                          	// j usually indicates decomposition scale
    int k;                                                          	// k is another variable used to denote spatial index

    //------- Declare collocation points -------------------------------//
    CollocationPoint** collPnt = new CollocationPoint*[J+1];            // J+1 rows
    for (int j=0;j<=J;j++) {                                            // 
        int N = jPnts(j);                                               // number of points at level j
        collPnt[j] = new CollocationPoint[N];                           // create objects for each of the points
    }                                                                   //

    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                         	    // 
        int N = ( jPnts(j) - 1 ) / 2;                                   //
        for (k=-N;k<=N;k++) {                                          	//
            collPnt[j][k+N].x = 2. * pow(2.,-(j+shift)) * k;            // x-locations of each collocation point
        }                                                       	    //
    }                                                           	    //

    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(collPnt);                                                 // decompose integrand into scaling and detail coefficients 

    //------- Remove coefficients below the threshold ------------------//
    thresholding(collPnt,threshold);                                    // knock out small wavelet coefficients

    //------- Compute integrals of wavelets at all levels --------------//
    wavelet_integrals(collPnt);                                         // compute integrals of the wavelets

    //------- Reconstruct function using wavelets ----------------------//    
    cout << integrate(collPnt)  << endl;                                               // compute integral of input function using wavelets

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;                                                   // delete collocation points from memory
    return 0;                                                           // 
}
