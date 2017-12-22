#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"

void time_integrate(CollocationPoint** collPnt,double h,string &equation, double c,double alpha) {

    //------- Advance solution in time ---------------------------------------------//
    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        collPnt[j][0].u = 1.;
        for (int i=1;i<N-1;i++) {
            if ( collPnt[j][i].isMask == true ) {
               collPnt[j][i].u = ( - ( collPnt[j][i].u + c ) * collPnt[j][i].ux + alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
            }
        }   
        collPnt[j][N-1].u = -1.;
    }
    return;
}
/*
        if ( equation.compare(0,7,"burgers") == 0 ) {
            collPnt[j][i].u = ( - collPnt[j][i].u * collPnt[j][i].ux + alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
        }

        if ( equation.compare(0,16,"modified_burgers") == 0 ) {
            collPnt[j][i].u = ( - collPnt[j][i].u * collPnt[j][i].ux + alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
        } */
