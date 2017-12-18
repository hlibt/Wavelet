#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void reconstruction_check(CollocationPoint** collPnt) {

    //--------------------------------------------------------------------------//
    // Information: reconstruction_check.cpp is called after both the 
    //              thresholding.cpp and adjacent_zone.cpp routines are ran.
    //              Once both routines have produced a mask of points which are
    //              currently significant (thresholding.cpp) and also points which
    //              may become significant at the next timestep (adjacent_zone.cpp),
    //              the final step to a complete mask is to ensure that the
    //              interpolation stencil used to compute each active wavelet
    //              coefficient is also included in the mask. These points may 
    //              include points corresponding to both inactive wavelets and 
    //              scaling coefficients alike. This step is known as the 
    //              perfect reconstruction check.
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

    //------- Loop through all points at levels j=2 through j=J ------------//
    for (int j=2;j<=J;j++) {                                                //
        int N = jPnts(j);                                                   //
        for (int k=0;k<N;k++) {                                             //
            if ( collPnt[j][k].isMask == true &&                            //
                    collPnt[j][k].isOdd == true) {                          //

                //------- Place interpolation stencil into mask --- --------//
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
                        collPnt[j-1][l].isMask = true;                      // place each stencil point in the mask
                }                                                           //

            }                                                               // 
        }                                                                   //
    }                                                                       //
}
