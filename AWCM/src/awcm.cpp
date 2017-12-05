#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include "CollocationPoint.hpp"
#include "matrix_solvers/BiCGSTAB.h"
#include "conditions/conditions.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
#include "global.hpp"
using namespace std;
void output(char* input, CollocationPoint** collPnt);

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
    J = 4;                                                              // number of scales in the system
    interpPnts = 2;                                                     // half the number of points used for interpolation (2*interpPnts + 1)
    double threshold = 0.;                                         	    // error tolerance for wavelet coefficients (determines accuracy of solution)
    int num_active = jPnts(J);                                          // declare variable to count number of active wavelets
    int i;                                                              // the usual counter variable for spatial index
    int j;                                                          	// j usually indicates decomposition scale
    int k;                                                          	// k is another variable used to denote spatial index

    //------- Physical parameters --------------------------------------//
    double advec_vel = 0.2;                                             // advection velocity

    //------- Define timestep size -------------------------------------//
    int num_steps = 1000;                                         	    // number of timesteps     
    double ti = 0.;                                               	    // initial simulation time  
    double tf = 0.01;                                              	    // final simulation time     
    double dt = ( tf - ti ) / num_steps;                           	    // timestep size (determined by choice of interval and number of steps)

    //------- Declare collocation points -------------------------------//
    CollocationPoint** collPnt = new CollocationPoint*[J+1];            // Create J+1 pointers to arrays
    for (int j=0;j<=J;j++) {                                            // 
        int N = jPnts(j);                                               // number of points at level j
        collPnt[j] = new CollocationPoint[N];                           // create array for all points at level j
    }                                                                   //
                        
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                         	    // 
        int N = ( jPnts(j) - 1 ) / 2;                                   //
        for (k=-N;k<=N;k++) {                                          	//
            collPnt[j][k+N].x = 2. * pow(2.,-(j+shift)) * k;            // x-locations of each collocation point
        }                                                       	    //
    }                                                           	    //

    //------- Sample the the initial condition -------------------------//
    for (k=0;k<jPnts(J);k++) {                                          //
        collPnt[J][k].scaling_coeff = init_condition(collPnt[J][k].x);  // sample the initial condition on the finest scale 
    }                                                                   //
    
    //------- Iterate in time ------------------------------------------//
    for (int t=0;t<num_steps;t++) {

        //------- Set scaling coefficients to new solution -------------//
        if ( t > 0 ) {
            for (i=0;i<jPnts(J);i++) {
                collPnt[J][i].scaling_coeff = collPnt[J][i].u;
            }
        }

        //------- Perform forward wavelet transform --------------------//
        fwd_trans(collPnt);                                             // compute all scaling and detail coefficients

        //------- Remove coefficients below the threshold --------------//
        thresholding(collPnt,threshold);                                // knock out small wavelet coefficients

        //------- Extend mask to include adjacent zone -----------------//

        //------- Compute spatial derivatives --------------------------//
        compute_derivatives(collPnt);                                   // compute the spatial derivatives

        //------- Reconstruct function using wavelets ------------------//    
        reconstruction(collPnt);                                        // build the solution using wavelets

        //------- Output solution at each timestep ---------------------//
        char u_out[45];
        sprintf(u_out,"output/_soln_files/u%d.dat",t);
        output(u_out,collPnt);

        //------- Advance in time --------------------------------------//
        time_integrate(collPnt,dt,advec_vel);                           // advance in time

    }                                                                   // end of time integration

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;
    return 0; 
}

void output(char* input, CollocationPoint** collPnt) {
    ofstream output;                            	
    output.open(input);               
    int k = 0;
    int cntr = 0;    
    bool* active = new bool[jPnts(J)];
    for (int i=0;i<jPnts(0);i++) { 
        k = indexShift(J,0,i);
        active[k] = true;
        cntr++;
    }
    for (int j=1;j<=J;j++) {
        for (int i=0;i<jPnts(j);i++) {
            if ( collPnt[j][i].isMask == true ) {
                k = indexShift(J,j,i);
                active[k] = true;
                cntr++;
            }
        }
    }
    double* x_compressed = new double[cntr+1];
    double* u_compressed = new double[cntr+1];
    k = 0;
    for (int i=0;i<jPnts(J);i++) {
        if ( active[i] == true ) {
            x_compressed[k] = collPnt[J][i].x;
            u_compressed[k] = collPnt[J][i].u;
            k++;
        }
    }
    for (int i=0;i<=cntr;i++) {  
 	    output << std::fixed << std::setprecision(16) << x_compressed[i] << " " << u_compressed[i] << endl;
    }
    output.close(); 
    delete[] active;
    delete[] x_compressed;
    delete[] u_compressed;
}
