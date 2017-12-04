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
using namespace std;
void output(char* input, double* U, CollocationPoint** collPnt);
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

//------- inline function declarations ---------------------------------//              
int inline jPnts(int j) { return pow(2,j+shift) + 1; }                  // the number of grid points at level j
                                                                        //
int inline indexShift(int jstar, int j, int k) {                        // represents index k at level j, at the desired level jstar
        int gridMultplr = pow(2,jstar-j);                               // multiplier parameter to get k to level jstar
        return k * gridMultplr;                                         // output what k would be at finest level jstar
}                                                                       //
                                                                        //
double inline init_condition(double x) {                                //
    return sin(2.*x);                                                   //
}                                                                       //


int main(void) {

    //------- Grid and tolerance parameters ----------------------------//
    shift = 2;                                                          // increases number of points on level j=0 (global variable)
    J = 12;                                                             // number of scales in the system
    interpPnts = 2;                                                     // half the number of points used for interpolation (2*interpPnts + 1)
    double threshold = pow(10.,-6);                                	    // error tolerance for wavelet coefficients (determines accuracy of solution)
    int i;                                                              // the usual counter variable for spatial index
    int j;                                                          	// j usually indicates decomposition scale
    int k;                                                          	// k is another variable used to denote spatial index

    //------- Define timestep size -------------------------------------//
    int num_steps = 40;                                            	    // number of timesteps     
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

    //------- sample the the initial condition -------------------------//
    for (k=0;k<jPnts(J);k++) {                                          //
        collPnt[J][k].scaling_coeff = init_condition(collPnt[J][k].x);  // sample the initial condition on the finest scale 
    }                                                                   //
    
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(collPnt);                                                 // compute all scaling and detail coefficients

    //------- Remove coefficients below the threshold ------------------//
    num_active = thresholding(collPnt,threshold);                       // knock out small wavelet coefficients

    //------- Extend mask ----------------------------------------------//

    //------- Compute spatial derivatives ------------------------------//
/*    for (j=0;j<J;j++) {                                                 // 
        int N=jPnts(j);                                                 // number of points at current level
        int gridMultplr=pow(2,J-j);                                     // constant needed to get to same point at higher level
        if (mask[j+1][1]==false && activPnt[0]==false) {                // left boundary
            double xEval=x[j][0];                                       //
            ux[0]=lagrInterpD1(xEval,x[j],c[j],0,interpPnts,N);         // compute first derivative at left boundary
            uxx[0]=lagrInterpD2(xEval,x[j],c[j],0,interpPnts,N);        // compute second derivative
            activPnt[0]=true;                                           //
        }                                                               //
        for (k=0;k<N;k++) {                                             //
            if (mask[j+1][2*k+1]==false && mask[j+1][2*k-1]==false      // check if function is well approximated
                    && mask[j][k]==true                                 // check if point is in the mask
                    && activPnt[gridMultplr*k]==false) {                // ensure derivative not already computed
                double xEval=x[j][k];                                   // evaluation point
                ux[gridMultplr*k]=lagrInterpD1(xEval,x[j],c[j],k,       // compute derivative from lagrange polynomial
                                    interpPnts,N);                      //
                uxx[gridMultplr*k]=lagrInterpD2(xEval,x[j],c[j],k,      // compute second derivative from lagrange polynomial
                                    interpPnts,N);                      // 
                activPnt[gridMultplr*k]=true;                           // represent this point at solution time
            }                                                           //
        }                                                               //
        if (mask[j+1][N-2]==false                                       // check if function well approximated
                && activPnt[gridMultplr*(N-1)]==false) {                // right boundary
            double xEval=x[j][N-1];                                     // evaluation point
            ux[gridMultplr*(N-1)]=lagrInterpD1(xEval,x[j],c[j],         // compute first derivative at right boundary
                        N-1,interpPnts,N);                              //
            uxx[gridMultplr*(N-1)]=lagrInterpD2(xEval,x[j],c[j],        // compute second derivative
                        N-1,interpPnts,N);                              //
            activPnt[gridMultplr*(N-1)]=true;                           //
        }                                                               //
    }                                                                   // */

    //------- Reconstruct function using wavelets ----------------------//    
    double* u = new double[num_active];                                 // vector containing the approximated function
    reconstruction(u,collPnt,num_active);                               // build the solution using wavelets
    
    //------- Output solution to file ----------------------------------//
    char u_Out[45];
    sprintf(u_Out,"output/_soln_files/u.dat");
    output(u_Out,u,collPnt);
//    coeffs_out(x,mask,J);

    //------- Setup rhs of pde -----------------------------------------//
/*    double* rhs=new double[jPnts(J)];                                   //
    double alp=0.00001;
    for (i=0;i<jPnts(J);i++) 
            rhs[i]=-u[i]*ux[i]+alp*uxx[i];
    } */

    //------- Advance in time ------------------------------------------//
//    RK4(u,rhs,activPnt,dt,jPnts(J));

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;
    return 0; 
}

void output(char* input, double* U, CollocationPoint** collPnt) {
    ofstream output;                            	
    output.open(input);                            	 
    for (int i=0;i<jPnts(J);i++) {  
 	    output << std::fixed << std::setprecision(16) << collPnt[J][i].x << " " << U[i] << endl;
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
