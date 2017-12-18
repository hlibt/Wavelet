#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"

void time_integrate(CollocationPoint** collPnt,double h,double c,double alpha) {

    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        collPnt[j][0].u = 0.;
        for (int i=1;i<N-1;i++) {
            if ( collPnt[j][i].isMask == true ) {
               collPnt[j][i].u = ( - collPnt[j][i].u * collPnt[j][i].ux + alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
            }
        }   
        collPnt[j][N-1].u = 0.;
    }
    return;
}
