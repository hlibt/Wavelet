#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
#include "../interpolation/interpolation.hpp"

void compute_derivatives(CollocationPoint** collPnt) { 

    //--------------------------------------------------------------------------------------//
    // Information: compute_derivative.cpp performs the forward wavelet                  
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.ux              - the derivative of 'u' to be computed 
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------------------//     

    //------- Compute spatial derivatives using finite differences -----//
    for (int i=0;i<jPnts(0);i++) {                                      // compute derivatives for all points associated with scaling functions
        double xeval = collPnt[0][i].x;                                 //
        collPnt[0][i].ux = lagrInterpD1(xeval,collPnt,J,i,jPnts(0));    // compute first derivatives on level j=0
    }                                                                   //
    for (int j=1;j<=J;j++) {                                            // 
        int N = jPnts(j);                                               // number of points at level j
        for (int i=0;i<N;i++) {                                         //
            if ( collPnt[j][i].isMask == true ) {                       // check if point is in the mask
                double xeval = collPnt[j][i].x;                         // evaluation point
                collPnt[j][i].ux = lagrInterpD1(xeval,collPnt,          // compute derivative from lagrange polynomial
                                        j-1,(i-1)/2,jPnts(j-1));        //
            }                                                           //
        }                                                               //
    }                                                                   //

    return;                                                 
}
