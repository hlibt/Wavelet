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
    for (int j=0;j<=J;j++) {                                                        // search all levels for points in the mask
        int N = jPnts(j);                                                           // number of points at level j
        for (int i=0;i<N;i++) {                                                     // loop through points at level j
 
            //------- Advance point if it is in mask -------------------------------//
            if ( collPnt[j][i].isMask == true ) {                                   // determine if point is in mask
                RK2(collPnt,j,i,h,equation);                                        // call time integration scheme
            }                                                                       //

            //------- Compute boundary values --------------------------------------//
            if ( boundary_type.compare(0,9,"derichlet") == 0 ) {
                collPnt[j][0].u = left_bc;
                collPnt[j][N-1].u = right_bc;
            }

        }   
    }
    return;
}
