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
    //              program does not ensure that the scaling coefficients used
    //              to compute each detail coefficient are kept, as this step
    //              is done in the reconstruction_check.cpp program.
    //              This collection of points needed to compute a detail  
    //              coefficient is called the interpolation stencil:
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
    
    //------- Decide fate of detail coefficients ---------------------------//
    for (int j=1;j<=J;j++) {                                                // loop through all other levels
        int N = jPnts(j);                                                   // number of points at level j
        for (int k=0;k<N;k++) {                                             // loop through points at level j

            //------- If coefficient less than threshold -------------------//
            if ( collPnt[j][k].isOdd == true &&                             // if point associated with wavelet but ...
                        abs(collPnt[j][k].detail_coeff) < epsilon ) {       // threshold not met, turn off the wavelet

                //------- Turn off wavelet ---------------------------------//
                collPnt[j][k].isMask = false;                               // turn off corresponding wavelet by removing point from mask
                collPnt[j][k].isBuffer = false;                             // remove possible previous affiliation with buffer zone

                //------- Remove neighbors from mask -----------------------//
                collPnt[j][k-1].isMask = false;                             // need to turn off left neighbor (scaling point)
                collPnt[j][k-1].isBuffer = false;                           // ensure removal from buffer  
                collPnt[j][k+1].isMask = false;                             // right neighbor as well (scaling point also)   
                collPnt[j][k+1].isBuffer = false;                           // ensure removal from buffer 
            }                                                               // 

            //------- If coefficient greater or equal to threshold ---------//
            if ( collPnt[j][k].isOdd == true &&                             // check if point corresponds to a wavelet
                 abs(collPnt[j][k].detail_coeff) >= epsilon ) {             // determine if detail coefficient is large enough

                //------- Place the point in the mask ----------------------//
                collPnt[j][k].isMask = true;                                // put active wavelet into the mask
                collPnt[j][k].isBuffer = false;                             // remove possible previous affiliation with adjacent zone

            }                                                               //

        }                                                                   //
    }                                                                       // 

    return;
}
