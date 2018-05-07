#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void decompress(CollocationPoint** collPnt, double* funcpts, int* spatial_identifier, int* level_identifier, int nactive) {

    //--------------------------------------------------------------------------//
    // Information: This program converts the compressed grid 
    //              back to the dyadic grid for the wavelet transform.
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
    
    //------- decompress the grid ------------------------------------------//
    for (int i=0;i<nactive;i++) {                                           //
        int j = level_identifier[i];                                        // dyadic grid level
        int k = spatial_identifier[i];                                      //
        collPnt[j][k].u = funcpts[i];                                       // updated solution on dyadic grid
    }                                                                       //

    //------- complete update ----------------------------------------------//
    for (int j=0;j<J;j++) {                                                 //
        int N = jPnts(j);                                                   //
        for (int k=0;k<N;k++) {                                             //
            if ( collPnt[j][k].isMask == true &&                            //
                    ( collPnt[j][k].isOdd == true || j==0 ) ) {             //
                collPnt[j+1][2*k].u = collPnt[j][k].u;                      //
            }                                                               // 
        }                                                                   //
    }                                                                       //
    return;
}
