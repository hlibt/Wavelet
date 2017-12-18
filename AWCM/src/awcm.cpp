#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "CollocationPoint.hpp"
#include "matrix_solvers/BiCGSTAB.h"
#include "conditions/conditions.hpp"
#include "transform/transform.hpp"
#include "grid_adaptation/grid_adaptation.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
#include "global.hpp"
#include "output/output.hpp"
using namespace std;
    
void control(int &max_scale, int &shift, double &threshold, int &interp_points, int &num_timesteps, double &tf, 
                double &advec_vel, double &diffusivity, string &buffer_type, int &buffer_width,
               int &buffer_height,  bool &ifwrite);

    //------------------------------------------------------------------------------//
    //                                                                              //
    //     ADAPTIVE WAVELET COLLOCATION METHOD USING 2ND GENERATION WAVELETS TO     //
    //     SOLVE THE 1D VISCOUS BURGERS EQUATION ON AN ADAPTIVE-DYADIC GRID         //
    //                                                                              //
    //     AUTHOR: BRANDON GUSTO                                                    //
    //     REFERENCES: O. Vasilyev & C. Bowman, JCP 165, 2000                       //
    //     DATE CREATED: Aug 01, 2017                                               //
    //                                                                              //
    //------------------------------------------------------------------------------//

int shift;                                                              // starting grid level
int J;                                                                  // maximum wavelet level    
int interpPnts;                                                         // half the number of interpolation points in each stencil

int main(void) {

    //------- Grid and tolerance parameters ----------------------------//
    double threshold;                                    	            // error tolerance for wavelet coefficients
    int buffer_width;                                                   // width of the buffer layer ( in terms of adjacent wavelets )
    int buffer_height;                                                  // height of the buffer layer ( in terms of wavelet levels )
    int i, j, k;                                                        // j is the wavelet level, i or k denote spatial index

    //------- Physical parameters --------------------------------------//
    double advec_vel;                                                   // advection velocity
    double diffusivity;                                                 // coefficient of diffusivity

    //------- Define timestep size -------------------------------------//
    int num_steps;                                               	    // number of timesteps     
    double tf;                                                    	    // final simulation time     

    //------- Numerical scheme choices ---------------------------------//
    string time_integrator;                                             // choice of time-stepping scheme
    string buffer_type;                                                 // defines how to define the adjacent zone layer of wavelets

    //------- Input/output preferences ---------------------------------//
    bool ifwrite;                                                       // write data to file or not

    //------- Input simulation parameters from control file ------------//
    control(J, shift, threshold, interpPnts, num_steps, tf,             // read all simulation input variables from file
               advec_vel, diffusivity, buffer_type, buffer_width,       //
               buffer_height, ifwrite);                                 //

    //------- Compute simulation timestep ------------------------------//
    double dt = tf / num_steps;                                    	    //

    //------- Allocate space for entire grid ---------------------------//
    CollocationPoint** collPnt = new CollocationPoint*[J+1];            // Create J+1 pointers to arrays
    for (int j=0;j<=J;j++) {                                            // 
        int N = jPnts(j);                                               // number of points at level j
        collPnt[j] = new CollocationPoint[N];                           // create array for all points at level j
    }                                                                   //

    //------- Populate dyadic grid -------------------------------------//
    seed_grid(collPnt);                                                 // sample initial condition on the grid

    //------- Iterate in time ------------------------------------------//
    for (int t=0;t<num_steps;t++) {                                     // step through time

        //------- Set scaling coefficients to new solution -------------//
        fwd_trans(collPnt);                                             // compute forward wavelet transform on the adaptive grid    

        //------- Remove coefficients below the threshold --------------//
        thresholding(collPnt,threshold);                                // construct an initial mask of wavelets whose coefficients are above the threshold

        //------- Extend mask to include adjacent zone -----------------//
        adjacent_zone(collPnt,buffer_width,buffer_height);              // create an adjacent zone of wavelets which may become significant after dt time

        //------- Perform the perfect reconstruction check -------------//
        reconstruction_check(collPnt);                                  // ensure that all detail points can be reconstructed at the next timestep

        //------- Reconstruct function using wavelets ------------------//    
        compute_field(collPnt);                                         // compute the field where necessary using wavelet basis

        //------- Output solution at each timestep ---------------------//
        if ( ifwrite == 1 ) write_field_to_file(collPnt,t);             // output solution to file at current timestep

        //------- Advance in time --------------------------------------//
        time_integrate(collPnt,dt,advec_vel,diffusivity);               // advance the solution forward in time

    }                                                                   // end of time integration

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;
    return 0; 
}

