#include <cmath>
#include <string>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
#include "phys_driver.hpp"

void RK2(double* gridpts, double* funcpts, int nactive, double h, double alp, double c, string equation) {

    //--------------------------------------------------------------------------//
    // Information: RK2.cpp updates grid points based on the second-order
    //              Runge-Kutta method. Each call to this program updates a 
    //              single grid point. The grid level j, and spatial index i
    //              of the grid point are input to the program.
    //
    // Input: 
    //              gridpts                 - the non-uniform grid of spatial points
    //              funcpts                 - functional value at non-uniform points 
    //              i                       - spatial index
    //              h                       - the size of the timestep
    //              equation                - a string to inform solver of pde
    // Output:
    //              collPnt.u               - the solution u_ji at the updated time
    //--------------------------------------------------------------------------//     
    
    for (int i=0;i<nactive;i++) {

        //------- Declare variables ------------------------        //
        double k1, k2;                                              // slopes of function at interval times

        //------- Compute RHS at each time interval --------        //
        k1 = rhs(gridpts,funcpts,nactive,i,0.,alp,c,equation);      // halfstep
        k2 = rhs(gridpts,funcpts,nactive,i,h*k1/2.,alp,c,equation); // modified slope

        //-------- Update the solution variable ------------        //
        funcpts[i] += h*k2;                                         // updated solution

    }

    return;
}
