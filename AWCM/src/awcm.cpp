#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "CollocationPoint.hpp"
#include "matrix_solvers/BiCGSTAB.h"
#include "conditions/conditions.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
using namespace std;
void output(char* input,double* X,double* U,bool* ACTIVE,int Jmax);
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
    int numSteps = 40;                                            	    // number of timesteps     
    double ti = 0.;                                               	    // initial simulation time  
    double tf = 1.;                                               	    // final simulation time     
    double dt = ( tf - ti ) / numSteps;                           	    // timestep size (determined by choice of interval and number of steps)

    //------- Class declarations ---------------------------------------//
    initial_condition IC;                                               // class containing initial conditions

    //------- Declare collocation points -------------------------------//
    int N = jPnts(0);                                                   // the number of points at the coarsest level j=0 (scaling points)
    ScalingPoint* scalPnt = new ScalingPoint[N];                        // declare points corresponding to scaling wavelets as objects
    DetailPoint** detPnt = new DetailPoint*[J+1];                       // J+1 rows
    for (int j=0;j<J;j++) {                                             // 
        N = jPnts(j) - 1;                                               // number of points associated with wavelets at level j
        detPnt[j] = new DetailPoint[N];                                 // create "DetailPoint" objects for each of the points associated with wavelets
    }                                                                   //
                        
    //------- Populate dyadic grid -------------------------------------//
    N = jPnts(0) - 1;                                                   // number of points at the coarsest level (-1)
    for (k=-N/2;k<=N/2;k++) {                                           //
        scalPnt[k+N].x = 2.*pow(2.,-(j+shift)) * k;                     // the x-values of scaling points on [-1,1]
    }                                                                   //
    for (j=0;j<J;j++) {                                         	    // 
        N = jPnts(j) - 1;                                               //
        for (k=-N;k<=N;k++) {                                          	//
            if ( k % 2 == 1 ) {                                         // determine if at odd point (wavelet location)
                detPnt[j][(k+N-1)/2].x = 2. * pow(2.,-(j+shift)) * k;   // x-locations of point corresponding to wavelet
            }                                                           //
        }                                                       	    //
    }                                                           	    //

    //------- sample the the initial condition -------------------------//
    for (k=0;k<jPnts(J);k++) {                                          //
        collPnt[J][k].scaling_coeff = IC.f( collPnt[J][k].x );          // sample the initial condition on the finest scale available
    }                                                                   //

    //------- Perform forward wavelet transform ------------------------//
//    fwd_trans(collPnt,J,interpPnts);                                    // decompose 'signal' into scaling and detail coefficients 
                                                                        // ... lifting not implemented yet

    //------- Remove coefficients below the threshold ------------------//
//    thresholding(collPnt,threshold,J);                                  // knock out small wavelet coefficients

    //------- Extend mask ----------------------------------------------//
/*    for (j=1;j<=J;j++) {                                                // 
        int jstar=j-1;                                                  //
        int gridMultplr=pow(2,j-jstar);                                 //
        for (k=0;k<jPnts(jstar);k++) {                                  //
            if (mask[jstar][k]==true) {                                 //
                mask[j][k*gridMultplr]=true;                            //
            }                                                           //
        }                                                               //
    }                                                                   //
    for (j=J-1;j>0;j--) {                                               // 
        for (k=0;k<jPnts(j);k++) {                                      //
            if (k%2==1 && mask[j][k]==true) {                           //
                int leftPnt=-interpPnts+1;  
                int rightPnt=interpPnts;    
                while ( leftPnt < 0 ) {
                    leftPnt++;
                    rightPnt++;
                }
                while ( rightPnt > (jPnts(j-1)-1) ) {
                    leftPnt--;
                    rightPnt--;
                }
                for (int l=leftPnt;l<=rightPnt;l++) {                   //
                    mask[j-1][(k-1)/2+l]=true;                          //
                }                                                       //
            }                                                           //
        }                                                               //
    }                                                                   //
*/
/*    //------- Calculate spatial derivatives ----------------------------//
    for (j=0;j<J;j++) {                                                 // 
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
//    reconstruction(x,u,c,d,mask,J,interpPnts);                          // build the solution using wavelets

    //------- Output solution to file ----------------------------------//
/*    char u_Out[45];
    char ux_Out[45];
    char uxx_Out[45];
    sprintf(u_Out, "output/_soln_files/u_t%d.dat",t);
    sprintf(ux_Out, "output/_soln_files/ux_t%d.dat",t);
    sprintf(uxx_Out, "output/_soln_files/uxx_t%d.dat",t);
    output(u_Out,x[J],u,activPnt,J);
    output(ux_Out,x[J],ux,activPnt,J);
    output(uxx_Out,x[J],uxx,activPnt,J);
    coeffs_out(x,mask,J); */

    //------- Setup rhs of pde -----------------------------------------//
/*    double* rhs=new double[jPnts(J)];                                   //
    double alp=0.00001;
    for (i=0;i<jPnts(J);i++) {
            rhs[i]=-u[i]*ux[i]+alp*uxx[i];
    } */

    //------- Advance in time ------------------------------------------//
//    RK4(u,rhs,activPnt,dt,jPnts(J));

    //------- Cleanup --------------------------------------------------//
    delete[] collPnt;
    return 0; 
}

void output(char* input,double* X,double* U,bool* ACTIVE,int Jmax) {
    ofstream output;                            	
    output.open(input);                            	 
    for (int i=0;i<jPnts(Jmax);i++) {  
        if (ACTIVE[i]==true) {
 			output << std::fixed << std::setprecision(16) << X[i] << " " << U[i] << endl;
        }
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
