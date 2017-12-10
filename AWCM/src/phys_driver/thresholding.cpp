#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void thresholding(CollocationPoint** collPnt, double epsilon) {

    //--------------------------------------------------------------------------//
    // Information: thresholding.cpp constructs an initial mask of points
    //              necessary for resolving the physics. These points will be
    //              integrated in time. This code however will not complete 
    //              the mask. A buffer zone of points must still be added to
    //              the mask so that physical features which become relavant
    //              after timestepping are not missed or underresolved. This 
    //              program also ensures that the scaling coefficients used
    //              to compute each detail coefficient which is kept, are also
    //              kept. This collection of points corresponding to 
    //              scaling coefficients is called the interpolation stencil:
    //
    //              detail coefficient at level j                   --->      o
    //              interpolation stencil on j-1 for interpPnts=2   ---> o  o   o  o
    //
    //           
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              epsilon                 - threshold parameter for detail coeffs.
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------//     
    
    //------- Place all scaling points at level j=0 into the mask ----------//
    for (int i=0;i<jPnts(0);i++) {                                          // loop through points at leve j=0
        collPnt[0][i].isMask = true;                                        // put these points into the mask
    }                                                                       //

    //------- Evaluate all detail coefficients -----------------------------//
    for (int j=1;j<=J;j++) {                                                // loop through all other levels
        int N = jPnts(j);                                                   // number of points at level j
        for (int k=0;k<N;k++) {                                             // loop through points at level j
            if ( collPnt[j][k].isOdd == true &&                             // check if point corresponds to a wavelet
                 abs(collPnt[j][k].detail_coeff) >= epsilon ) {             // determine if detail coefficient is large enough

                //------- Place the point in the mask ----------------------//
                collPnt[j][k].isMask = true;                                // put active wavelet into the mask

                //------- Extend mask to its interpolation stencil ---------//
                int leftPnt = -interpPnts + 1 + (k-1)/2;                    // determine the leftmost point in its stencil
                int rightPnt = interpPnts + (k-1)/2;                        // determine the rightmost point in the stencil
                while ( leftPnt < 0 ) {                                     // this loop checks to ensure the stencil consists of points in the domain
                    leftPnt++;                                              // adjust stencil one point to the right if necessary
                    rightPnt++;                                             // 
                }                                                           //
                while ( rightPnt > (jPnts(j-1)-1) ) {                       // adjust stencil one point to the left if necessary
                    leftPnt--;                                              //
                    rightPnt--;                                             //
                }                                                           //
                for (int l=leftPnt;l<=rightPnt;l++) {                       // loop through points in the stencil
                    collPnt[j-1][l].isMask = true;                          // place each point in the interpolation stencil in the mask ...
                }                                                           // so that the detail coefficient can be reconstructed ...
            }                                                               // at the next timestep
            else if ( collPnt[j][k].isOdd &&                                // if point associated with wavelet but ...
                        abs(collPnt[j][k].detail_coeff) < epsilon ) {       // threshold not met, turn off the wavelet
                    collPnt[j][k].isMask = false;                           //
            }                                                               //
        }                                                                   //
    }                                                                       // a near-complete mask is constructed ( buffer layer still needed )

    return;                                                                 
}
