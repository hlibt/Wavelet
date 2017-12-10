#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "CollocationPoint.hpp"
#include "matrix_solvers/BiCGSTAB.h"
#include "conditions/conditions.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
#include "global.hpp"
using namespace std;
void output(CollocationPoint** collPnt,int current_timestep);

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

int shift;
int J;
int interpPnts;

int main(void) {

    //------- Grid and tolerance parameters ----------------------------//
    shift = 2;                                                          // increases number of points on level j=0 (global variable)
    J = 7;                                                              // number of scales in the system
    interpPnts = 2;                                                     // half the number of points used for interpolation (2*interpPnts + 1)
    double threshold = pow(10.,-6.);                         	        // error tolerance for wavelet coefficients (determines accuracy of solution)
    int i, j, k;                                                        // j is the wavelet level, i or k denote spatial index

    //------- Physical parameters --------------------------------------//
    double advec_vel = 1. ;                                             // advection velocity

    //------- Define timestep size -------------------------------------//
    int num_steps = 10000;                                         	    // number of timesteps     
    double ti = 0.;                                               	    // initial simulation time  
    double tf = 0.5;                                              	    // final simulation time     
    double dt = ( tf - ti ) / num_steps;                           	    // timestep size (determined by choice of interval and number of steps)

    //------- Declare collocation points -------------------------------//
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
        thresholding(collPnt,threshold);                                // construct a semi-complete mask for the next timestep

        //------- Extend mask to include adjacent zone -----------------//
        extend_mask(collPnt,1);                                         // extend the mask to include all points necessary for computing detail coeffs

        //------- Compute spatial derivatives --------------------------//
//        compute_derivatives(collPnt);                                   // compute the spatial derivatives

        //------- Reconstruct function using wavelets ------------------//    
        reconstruction(collPnt);                                        // build the solution using wavelets

        //------- Output solution at each timestep ---------------------//
        output(collPnt,t);                                              // output solution to file at current timestep

        //------- Advance in time --------------------------------------//
        time_integrate(collPnt,dt,advec_vel);                           // advance the solution forward in time

    }                                                                   // end of time integration

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;
    return 0; 
}

void output(CollocationPoint** collPnt,int current_timestep) {

    char u_out[45];                                                     // declare space for character array
    sprintf(u_out,"output/_soln_files/u%d.dat",current_timestep);       // append filename to include timestep number
    ofstream output;                            	
    output.open(u_out);     
    int num_active = 0;    
    CollocationPoint* tmp = new CollocationPoint[jPnts(J)];
    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {  
            if ( collPnt[j][i].isMask == true ) {
                int k = indexShift(J,j,i);
                tmp[k] = collPnt[j][i];
                num_active++;
            }
        }
    }
    for (int i=0;i<jPnts(J);i++) {
        if ( tmp[i].isMask == true ) {
            output << std::fixed << std::setprecision(16) << tmp[i].x << " " << tmp[i].u << endl;    
        }
    }
    output.close(); 
    printf(" \n");
    printf("------------------------------------------ \n");
    printf("Timestep: %d \n",current_timestep);
    printf("Active points: %d out of %d \n",num_active,jPnts(J));
    printf("------------------------------------------ \n");   
    printf(" \n");
    delete[] tmp;
}
