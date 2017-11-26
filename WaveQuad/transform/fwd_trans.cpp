#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
using namespace std;

void fwd_trans(CollocationPoint** collPnt) {
    
    //--------------------------------------------------------------------------//
    // Information: fwd_trans performs the forward wavelet                  
    //              transformation to compute the scaling and detail
    //              coefficients for all levels 
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interpolation points
    //--------------------------------------------------------------------------//     
    
    for (int i=0;i<jPnts(J);i++) {                                              // loop through the points at level Jmax
        collPnt[J][i].scaling_coeff = integrand( collPnt[J][i].x );             // set the scaling coefficients at level J to the integrand
    };                                                                          //
    for (int j=J;j>0;j--) {                                                     // copy even coefficients from level j to j-1
        int n = jPnts(j);                                                       // 
        for (int i=0;i<n;i++) {                                                 // loop through points at level j
            if (i%2==0) {                                                       //
            collPnt[j-1][i/2].scaling_coeff = collPnt[j][i].scaling_coeff;      // even scaling coeff's the same, just copy them
            }                                                                   //
        }                                                                       //
    }                                                                           //
    for (int j=J;j>0;j--) {                                                     //
        int n = jPnts(j);                                                       //
        for (int i=0;i<n;i++) {                                                 //
            if ( i%2==1 ) {                                                     // detail coefficients only exist at odd points
            double xeval = collPnt[j][i].x;                                     // define the point for polynomials to be evaluated 
            collPnt[j][i].detail_coeff = .5*( collPnt[j][i].scaling_coeff       // compute detail coefficients
                        - lagrInterp(xeval,collPnt[j-1],(i-1)/2,jPnts(j-1)) );  //
            }                                                                   //
        }                                                                       //
    }                                                                           //
    return;                                                                     //
}  
