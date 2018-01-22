#include <cmath>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
#include "phys_driver.hpp"

void RK2(CollocationPoint** collPnt, int j, int i, double h, string equation) {

    //--------------------------------------------------------------------------//
    // Information: RK2.cpp updates grid points based on the second-order
    //              Runge-Kutta method. Each call to this program updates a 
    //              single grid point. The grid level j, and spatial index i
    //              of the grid point are input to the program.
    //
    // Input: 
    //              collPnt                 - the collocation point data structure
    //              j                       - grid level
    //              i                       - spatial index
    //              h                       - the size of the timestep
    // Output:
    //              collPnt.u               - the solution u_ji at the updated time
    //--------------------------------------------------------------------------//     
    
    //------- Declare variables ------------------------//
    double k1, k2;

    //------- Compute RHS at each time interval --------//
    k1 = rhs(collPnt,j,i,0.);
    k2 = rhs(collPnt,j,i,h*k1/2.);

    //-------- Update the solution variable ------------//
    collPnt[j][i].u = collPnt[j][i].u + h*k2;

    return
}
