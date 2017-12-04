#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
#include "../interpolation/interpolation.hpp"

void compute_derivatives( CollocationPoint** collPnt ) { 

    for (int i=0;i<jPnts(0);i++) {                                      // compute derivatives for all points associated with scaling functions
        double xeval = collPnt[0][i].x;                                 //
        int k = indexShift(J,0,i);                                      // constant needed to get to same point at higher level
        collPnt[J][k].ux = lagrInterpD1(xeval,collPnt,0,i,jPnts(0));    // compute first derivatives on level j=0
    }                                                                   //
    for (int j=1;j<=J;j++) {                                            // 
        int N = jPnts(j);                                               // number of points at level j
        for (int i=0;i<N;i++) {                                         //
            if ( collPnt[j][i].isMask == true ) {                       // check if point is in the mask
                double xeval = collPnt[j][i].x;                         // evaluation point
                int k = indexShift(J,j,i);                              // 
                collPnt[J][k].ux = lagrInterpD1(xeval,collPnt,          // compute derivative from lagrange polynomial
                                        j-1,(i-1)/2,jPnts(j-1));        //
            }                                                           //
        }                                                               //
    }                                                                   //

}
