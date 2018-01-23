#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void time_advance(CollocationPoint** collPnt,double h,string equation, double c,double alpha,
                        string boundary_type, double left_bc, double right_bc) {

    //------- Advance interior points in time --------------------------------------//
    for (int j=1;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {

            //------- Advance point if it is in mask -------------------------------//
            if ( collPnt[j][i].isMask == true && collPnt[j][i].isOdd == true ) {
                RK2(collPnt,j,i,h,equation);
            }

            //------- Compute boundary values --------------------------------------//
            if ( boundary_type.compare(0,9,"derichlet") == 0 ) {
                collPnt[j][0].u = left_bc;
                collPnt[j][N-1].u = right_bc;
            }

        }   
    }
    return;
}
