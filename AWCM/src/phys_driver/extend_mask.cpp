#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void extend_mask(CollocationPoint** collPnt,int buffer_width) {
    
    //--------------------------------------------------------------------------//
    // Information: extend_mask.cpp extends the existing mask to include
    //              a buffer zone, or a collection of points corresponding
    //              to wavelets which may become significant throughout the 
    //              next timestep. As well, the points in the interpolation
    //              stencil for those new wavelets are also included.
    //
    //           
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------//     
    
    //------- Extend mask to possibly significant detail coefficients ------//
    for (int j=1;j<=J;j++) {                                                // loop through all other levels
        int N = jPnts(j);                                                   // number of points at level j
        for (int k=0;k<N;k++) {                                             // loop through points at level j
            if ( collPnt[j][k].isOdd ==true &&                              // check if point corresponds to a wavelet
                    collPnt[j][k].isMask == true ) {                        // check if wavelet point is in the mask already

                //------- Check if buffer_width is still in the domain -----//
                if ( k + buffer_width*2 < N ) {                             //

                    //------- Place the point in the mask ------------------//
                    int knew = k + buffer_width * 2;                        // new index of the safety wavelet
                    collPnt[j][knew].isMask = true;                         // put new safety wavelet into the mask

                    //------- Extend mask to its interpolation stencil -----//
                    int leftPnt = -interpPnts + 1 + (knew-1)/2;             // determine the leftmost point in its stencil
                    int rightPnt = interpPnts + (knew-1)/2;                 // determine the rightmost point in the stencil
                    while ( leftPnt < 0 ) {                                 // this loop checks to ensure the stencil consists of points in the domain
                        leftPnt++;                                          // adjust stencil one point to the right if necessary
                        rightPnt++;                                         // 
                    }                                                       //
                    while ( rightPnt > (jPnts(j-1)-1) ) {                   // adjust stencil one point to the left if necessary
                        leftPnt--;                                          //
                        rightPnt--;                                         //
                    }                                                       //
                    for (int l=leftPnt;l<=rightPnt;l++) {                   // loop through points in the stencil
                        collPnt[j-1][l].isMask = true;                      // place each point in the interpolation stencil in the mask
                    }                                                       // 
                }                                                           //
            }                                                               // 
        }                                                                   //
    }                                                                       // a complete mask is now constructed

    return;
}
