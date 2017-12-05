#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"

void time_integrate(CollocationPoint** collPnt,double h,double c) {

/*    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {
            if ( collPnt[j][i].isMask == true ) {
                collPnt[j][i].u = ( -c * collPnt[j][i].ux ) * h + collPnt[j][i].u;
            }
        }   
    }
    collPnt[0][0] = collPnt[0][jPnts(0)-1]; */

    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {
            int k = indexShift(J,j,i);
            if ( collPnt[j][i].isMask == true ) {
                collPnt[J][k].u = ( -c * collPnt[J][k].ux ) * h + collPnt[J][k].u;
            }
        }   
    }
    collPnt[J][0].u = collPnt[J][jPnts(J)-1].u;
    return;
}
