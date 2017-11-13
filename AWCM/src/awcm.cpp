#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "matrix_solvers/BiCGSTAB.h"
#include "conditions/conditions.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
void output(char* input,double* X,double* U,bool* ACTIVE,int Jmax);
void coeffs_out(double** x,bool** mask,int Jmax);

    //------------------------------------------------------------------------------//
    //                                                                              //
    //     ADAPTIVE WAVELET COLLOCATION METHOD USING 2ND GENERATION WAVELETS TO     //
    //     SOLVE THE 1D VISCOUS BURGERS EQUATION ON AN ADAPTIVE-DYADIC GRID         //
    //                                                                              //
    //     AUTHOR: BRANDON GUSTO                                                    //
    //     DATE: Aug 01, 2017                                                       //
    //                                                                              //
    //------------------------------------------------------------------------------//

int shift;

// Inline functions:
int inline jPnts(int j) {return pow(2,j+shift)+1;}

int main(void) {

    //------- Grid and tolerance parameters ----------------------------//
    shift=2;                                                            // increases number of points of level j=0
    int J=12;                                                           // number of scales in the system
    int interpPnts=2;                                                   // half the number of points used for interpolation
    double threshold=pow(10.,-6);                                	    // error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j usually indicates decomposition scale
    int k;                                                          	// k is another variable for spatial index

    //------- Define timestep size -------------------------------------//
    int numsteps=40;                                            	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/numsteps;                                 	    // timestep size

    //------- Class declarations ---------------------------------------//
    initial_condition initCondition;                                    // class containing initial conditions

    //------- Declare collocation points -------------------------------//
    CollocationPoint** collPnt = new CollocationPoint*[J+1];            // collocation points as objects
    for (j=0;j<=J;j++) {                                                //
        int N = jPnts(j);                                               // number of points at level j
        collPnt[i] = new CollocationPoint[N];                           // 
    }                                                                   //
                        
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=jPnts(j-1)-1;                                             //
        for (k=-N;k<=N;k++) {                                    	    //
            collPnt[j][k+N].x = 2. * pow(2.,-(j+shift)) * k;            // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //

    //------- Initially set c's to initial condition -------------------//
    for (k=0;k<jPnts(J);k++) {                                          //
        collPnt[J][k].scaling_coeff = initCondition.f(x[J][k]);         // sample initial condition on finest scale 
    }                                                                   //

    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(collPnt,J,interpPnts);                                    // decompose 'signal' into scaling and detail coefficients

    //------- Remove coefficients below the threshold ------------------//
    thresholding(collPnt,threshold,J);                                  // knock out small wavelet coefficients

    //------- Extend mask ----------------------------------------------//
    for (k=0;k<jPnts(0);k++) mask[0][k]=true;                           // scaling coefficients at coarsest level included in mask
    for (j=1;j<=J;j++) {                                                // 
        int jstar=j-1;                                                  //
        int gridMultplr=pow(2,j-jstar);                                 //
        for (k=0;k<jPnts(jstar);k++) {                                  //
            if (mask[jstar][k]==true) {                                 //
                mask[j][k*gridMultplr]=true;                            //
            }                                                           //
        }                                                               //
    }                                                                   //
/*    for (j=J-1;j>0;j--) {                                               // 
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
    //------- Calculate spatial derivatives ----------------------------//
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
    }                                                                   //
    //------- Reconstruct function using wavelets ----------------------//    
    reconstruction(x,u,c,d,mask,J,interpPnts);                          // build the solution using wavelets

    //------- Output solution to file ----------------------------------//
    int t=1;
    char u_Out[45];
    char ux_Out[45];
    char uxx_Out[45];
    sprintf(u_Out, "output/_soln_files/u_t%d.dat",t);
    sprintf(ux_Out, "output/_soln_files/ux_t%d.dat",t);
    sprintf(uxx_Out, "output/_soln_files/uxx_t%d.dat",t);
    output(u_Out,x[J],u,activPnt,J);
    output(ux_Out,x[J],ux,activPnt,J);
    output(uxx_Out,x[J],uxx,activPnt,J);
    coeffs_out(x,mask,J);

    //------- Setup rhs of pde -----------------------------------------//
    double* rhs=new double[jPnts(J)];                                   //
    double alp=0.00001;
    for (i=0;i<jPnts(J);i++) {
            rhs[i]=-u[i]*ux[i]+alp*uxx[i];
    }

    //------- Advance in time ------------------------------------------//
//    RK4(u,rhs,activPnt,dt,jPnts(J));

    //------- Cleanup --------------------------------------------------//
    delete[] x;
    delete[] c;
    delete[] d;
    delete[] u;
    delete[] ux; 
    delete[] uxx;
    delete[] rhs;
    delete[] activPnt;
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
