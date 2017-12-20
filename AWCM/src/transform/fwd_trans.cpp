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
    // Information: fwd_trans.cpp performs the forward wavelet                  
    //              transform to compute the detail coefficients for
    //              points associated with currently active wavelets
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------//     
    
    //------- Loop through all levels, compute scaling coefficients ------------//
    for (int j=J;j>=0;j--) {                                                    // begin from finest level of resolution
        int n = jPnts(j);                                                       // number of points at level j
        for (int i=0;i<n;i++) {                                                 // loop through all points at level j
            if (collPnt[j][i].isMask == true) {                                 // check if point is in the mask from the last timestep
                collPnt[j][i].scaling_coeff = collPnt[j][i].u;                  // scaling coefficients are solution from last timestep
            }                                                                   //
        }                                                                       //
    }                                                                           //

    //------- Loop through all levels, compute detail coefficients -------------//
    for (int j=J;j>0;j--) {                                                     // begin from finest level of resolution
        int n = jPnts(j);                                                       // number of points at level j
        for (int i=0;i<n;i++) {                                                 // loop through all points at level j
            if (collPnt[j][i].isMask == true && collPnt[j][i].isOdd == true) {  // check if point associated with wavelet is active
                double xeval = collPnt[j][i].x;                                 // define the point for interpolating polynomial to be evaluated 
                collPnt[j][i].detail_coeff = .5 * (collPnt[j][i].scaling_coeff  // compute detail coefficients
                        - lagrInterp(xeval,collPnt,j-1,(i-1)/2,jPnts(j-1)));    //
            }                                                                   //
        }                                                                       //
    }                                                                           //

    return;                                                                     //
}  
