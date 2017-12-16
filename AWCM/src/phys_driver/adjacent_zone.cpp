#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void adjacent_zone(CollocationPoint** collPnt,int buffer_width,int buffer_height) {
    
    //--------------------------------------------------------------------------//
    // Information: adjacent_zone.cpp extends the existing mask to include
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
    //              collPnt.isMask          - boolean (in the current mask or not)
    //              buffer_width            - number of wavelets left and right to activate
    //              buffer_height           - number of levels above and below to activate
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------//     
    
    //------- Extend mask to possibly significant detail coefficients ------//
    for (int j=1;j<=J;j++) {                                                // loop through all other levels
        int N = jPnts(j);                                                   // number of points at level j
        for (int k=0;k<N;k++) {                                             // loop through points at level j
            if ( collPnt[j][k].isOdd == true && collPnt.isBuffer == false   // check if point corresponds to a wavelet which is not in safety zone
                    && collPnt[j][k].isMask == true ) {                     // check if wavelet point is in the mask already

                //---- Include points left and right of active wavelet -----//
                for (int l=2;l<=buffer_width*2;l=l+2) {                     // loop through 'buffer_width' number of points to the left
                    if ( k-l > 0 && collPnt[j][k-l].isMask == false ) {     // check if the point is within the grid and also not currently active 
                        collPnt[j][k-l].isMask = true;                      // put new safety wavelet into the mask
                        collPnt[j][k-l].isBuffer = true;                    // declare the safety wavelet as such, so it does not activate further wavelets 
                    }                                                       //
                    if ( k+l < N && collPnt[j][k+l].isMask == false ) {     // check if the point is within the grid and also not currently active 
                        collPnt[j][k+l].isMask = true;                      // put new safety wavelet into the mask
                        collPnt[j][k+l].isBuffer = true;                    // declare the safety wavelet as such, so it does not activate further wavelets 
                    }                                                       //
                }                                                           //

                //---- Include points above and below active wavelet -------//
                for (int jstar=1;jstar<=buffer_height;jstar++) {            // loop through 'buffer_height' number of wavelet levels
                    int kstar = indexShift(j+jstar,j,k);                    // compute index of same wavelet (at j)  at level 'jstar' of interest
                    if ( j + jstar <= J && kstar + 1 < jPnts(j+jstar) &&    // check if point level above and to the right is within the domain 
                            collPnt[j+jstar][kstar+1].isMask == false ) {   // check if this point is currently in the mask
                        collPnt[j+jstar][kstar+1].isMask = true;            // put point in the mask
                        collPnt[j+jstar][kstar+1].isBuffer = true;          // declare the safety wavelet as such, so it does not activate further wavelets
                    }                                                       //
                    if ( j + jstar <= J && kstar - 1 > 0 &&                 // check if point level above and to the left is within the domain 
                            collPnt[j+jstar][kstar-1].isMask == false ) {   // check if this point is currently in the mask
                        collPnt[j+jstar][kstar-1].isMask = true;            // put point in the mask
                        collPnt[j+jstar][kstar-1].isBuffer = true;          // declare the safety wavelet as such, so it does not activate further wavelets
                    }                                                       //
                    kstar = (k-1) / 2;                                      // wavelet point below and to the left of k
                    if ( j - jstar > 0 && kstar > 0 &&                      // check if point level above and to the left is within the domain 
                            collPnt[j-jstar][kstar].isOdd == true &&        // check if kstar is a wavelet point
                            collPnt[j-jstar][kstar].isMask == false ) {     // check if this point is currently in the mask
                        collPnt[j-jstar][kstar].isMask = true;              // put point in the mask
                        collPnt[j-jstar][kstar].isBuffer = true;            // declare the safety wavelet as such, so it does not activate further wavelets
                    }
                    kstar = (k+1) / 2;                                      // wavelet point below and to the left of k
                    if ( j - jstar > 0 && kstar < jPnts(j-jstar) &&         // check if point level above and to the right is within the domain 
                            collPnt[j-jstar][kstar].isOdd == true &&        // check if kstar is a wavelet point
                            collPnt[j-jstar][kstar].isMask == false ) {     // check if this point is currently in the mask
                        collPnt[j-jstar][kstar].isMask = true;              // put point in the mask
                        collPnt[j-jstar][kstar].isBuffer = true;            // declare the safety wavelet as such, so it does not activate further wavelets
                    }                                                       //
                }                                                           //
            }                                                               // 
        }                                                                   //
    }                                                                       // an adjacent zone is now constructed

    return;
}
