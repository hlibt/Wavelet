#include <iostream>
#include <cmath>
#include <math.h>
#include "../global.hpp"
#include "../CollocationPoint.hpp"
using namespace std;

double lagrInterp(double xeval, CollocationPoint** collPnt, int j, int i, int numPnts) {

    //--------------------------------------------------------------------------//
    // Information: This program "lagrInterp.cpp" computes the Lagrange interpolating
    //              polynomial, including the coefficients, and evaluates it at the 
    //              input variable x. This is used to construct wavelets, and in the 
    //              wavelet transform itself.
    //
    // Input: 
    //              xeval                    - the point that the polynomial is evaluated
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              J                       - maximum grid level (global)
    //              numPnts                 - the number of points at current level
    //--------------------------------------------------------------------------//     

    double sum = 0.;                                                            // initialize a summation variable
    int leftPnt = -interpPnts + 1 + i;                                          // first guess as to the leftmost point in the interp. stencil
    int rightPnt = interpPnts + i;                                              // same for the right point
    while ( leftPnt < 0 ) {                                                     // this loop checks to ensure the stencil consists of points in the domain
        leftPnt++;                                                              // adjust stencil one point to the right if necessary
        rightPnt++;                                                             // 
    }                                                                           //
    while ( rightPnt > (numPnts-1) ) {                                          // adjust stencil one point to the left if necessary
        leftPnt--;                                                              //
        rightPnt--;                                                             //
    }                                                                           //
    for (int l=leftPnt;l<=rightPnt;l++) {                                       // compute the weights
        double product = 1.;                                                    // 
        for (int k=leftPnt;k<=rightPnt;k++) {                                   //
            if (k==l) {                                                         // if l the same as k do nothing
            } else {                                                            // else compute the Lagrange coefficients
                product *= ( xeval - collPnt[j][k].x ) / ( collPnt[j][l].x      //
                                    - collPnt[j][k].x );                        //     
            }                                                                   //
        }                                                                       //
        sum += product * collPnt[j][l].scaling_coeff;                           // compute the polynomial
    }
    return sum;
}
