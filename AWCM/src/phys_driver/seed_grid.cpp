#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"

void seed_grid(CollocationPoint** collPnt) {

    //--------------------------------------------------------------------------//
    // Information: seed_grid.cpp initializes the collocation point objects
    //              and sets their locations on the grid, their solution to the
    //              initial condition, and whether the point is odd, implying 
    //              that it corresponds to a wavelet
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.u               - the solution variable
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------//     

    for (int j=0;j<=J;j++) {                                      	    // 
        int N = ( jPnts(j) - 1 ) / 2;                                   //
        for (int k=-N;k<=N;k++) {                                     	//
            double x = 2. * pow(2.,-(j+shift)) * k;                     // x-locations of each collocation point
            collPnt[j][k+N].x = x;                                      // populate the dyadic grid
            collPnt[j][k+N].u = init_condition(x);                      // sample the initial condition
            collPnt[j][k+N].isMask = true;                              // initially treat all points as in the mask
            if ( (k+N)%2==1 ) collPnt[j][k+N].isOdd = true;             // specifiy odd points (for wavelets/detail coefficients)
        }                                                       	    //
    }                                                           	    //
}
