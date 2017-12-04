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
void output(char* input, double* X, double* U, int N);
void coeffs_out(double** x,bool** mask,int Jmax);

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
    J = 10;                                                             // number of scales in the system
    interpPnts = 2;                                                     // half the number of points used for interpolation (2*interpPnts + 1)
    double threshold = pow(10.,-5);                                	    // error tolerance for wavelet coefficients (determines accuracy of solution)
    int num_active;                                                     // declare variable to count number of active wavelets
    int i;                                                              // the usual counter variable for spatial index
    int j;                                                          	// j usually indicates decomposition scale
    int k;                                                          	// k is another variable used to denote spatial index

    //------- Define timestep size -------------------------------------//
    int num_steps = 400;                                         	    // number of timesteps     
    double ti = 0.;                                               	    // initial simulation time  
    double tf = 1.;                                               	    // final simulation time     
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
            collPnt[j][k+N].x = 2. *  pow(2.,-(j+shift)) * k;           // x-locations of each collocation point
        }                                                       	    //
    }                                                           	    //

    //------- Sample the the initial condition -------------------------//
    for (k=0;k<jPnts(J);k++) {                                          //
        collPnt[J][k].scaling_coeff = init_condition(collPnt[J][k].x);  // sample the initial condition on the finest scale 
    }                                                                   //
    
    //------- Begin stepping in time -----------------------------------//
    for (int t=0;t<num_steps;t++) {

        //------- Sample the previous timestep -------------------------//

        //------- Perform forward wavelet transform --------------------//
        fwd_trans(collPnt);                                             // compute all scaling and detail coefficients

        //------- Remove coefficients below the threshold --------------//
        num_active = thresholding(collPnt,threshold);                   // knock out small wavelet coefficients

        //------- Extend mask ------------------------------------------//

        //------- Compute spatial derivatives --------------------------//
        compute_derivatives(collPnt);                                   // compute the spatial derivatives

        //------- Reconstruct function using wavelets ------------------//    
        reconstruction(collPnt);                                        // build the solution using wavelets

        //------- Compress the solution into active points -------------//
        bool* active = new bool[jPnts(J)];
        for (i=0;i<jPnts(0);i++) { 
            k = indexShift(J,0,i);
            active[k] = true;
        }
        for (j=1;j<=J;j++) {
            for (i=0;i<jPnts(j);i++) {
                if ( collPnt[j][i].isMask == true ) {
                    k = indexShift(J,j,i);
                    active[k] = true;
                }
            }
        }
        double* u_compressed = new double[ num_active ];
        double* ux_compressed = new double[ num_active ];
        double* x_compressed = new double[ num_active ];
        k = 0;
        for (i=0;i<jPnts(J);i++) {
            if ( active[i] == true ) {
                x_compressed[k] = collPnt[J][i].x;
                u_compressed[k] = collPnt[J][i].u;
                ux_compressed[k] = collPnt[J][i].ux;
                k++;
            }
        }
        delete[] active;
        
        //------- Output solution to file ----------------------------------//
        char u_out[45];
        char ux_out[45];
        sprintf(u_out,"output/_soln_files/u.dat");
        sprintf(ux_out,"output/_soln_files/ux.dat");
        output(u_out,x_compressed,u_compressed,num_active);
        output(ux_out,x_compressed,ux_compressed,num_active);

        //------- Setup rhs of pde -----------------------------------------//
        double* f = new double[num_active];                                 // right hand side vector
        double advec_vel = 0.2;                                             // advection velocity
        compute_rhs(f,ux_compressed,advec_vel,num_active);                  // compute the rhs of the pde

        //------- Advance in time ------------------------------------------//
        RK4(u_compressed, f, dt, num_active);                               // advance in time

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;
    delete[] x_compressed;
    delete[] u_compressed;
    delete[] ux_compressed;
    return 0; 
}

void output(char* input, double* X, double* U, int N) {
    ofstream output;                            	
    output.open(input);                            	 
    for (int i=0;i<N;i++) {  
 	    output << std::fixed << std::setprecision(16) << X[i] << " " << U[i] << endl;
    }
    output.close(); 
}

void coeffs_out(double** x,bool** mask,int Jmax) {
    for (int j=0;j<Jmax;j++) {
        int N=jPnts(j);   
        ofstream output;
        char fn[45];
        snprintf(fn,sizeof fn,"output/_coeff_files/coeff%d.dat",j);
        output.open(fn);
        if (j==0) {
            for (int t=0;t<N;t++) {
                output<<x[j][t]<<" "<<0<<endl;
            }
        } else {
            for (int t=0;t<N;t++) {
                if (mask[j][t]==true && t%2==1) {
                    output<<x[j][t]<<" "<<j<<endl;
                } else { 
                    output<<x[j][t]<<" "<<-10<<endl;
                }
            }
        }
        output.close();
    }
}
