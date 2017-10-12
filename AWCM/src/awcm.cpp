#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "matrix_solvers/BiCGSTAB.h"
#include "initial_conditions/initial_conditions.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
using namespace std;

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
    //------- General parameters ---------------------------------------//
    shift=2;                                                            // increases number of points of level j=0
    int J=10;                                                            // number of scales in the system
    int interpPnts=2;                                                   // half the number of points used for interpolation
    double threshold=pow(10.,-3);                               	    // error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    //------- Class definitions ----------------------------------------//
    initial_condition IC;                                           	// declare 'initial condition' class variable
    //------- Declare arrays -------------------------------------------//
    int N=jPnts(J);                                                     //
    double* u=new double[N];                      	                    // solution variable current timestep
    double* ux=new double[N];                                           // first derivative of the solution wrt x
    bool* activPnt=new bool[N];                                         // shows which points are active in grid
    double** x=new double*[J+1];			            		        // dyadic grid storage
    double** c=new double*[J+1];                                        // scaling function coefficients
    bool** mask=new bool*[J+1];                                         // mask denoting where detail coefficients are kept 
    double** d=new double*[J];                                          // detail function coefficients
    for (j=0;j<=J;j++) {						                        //
    	int N=jPnts(j); 			    	                            // number of points at level j
	    x[j]=new double[N];  						                    // dyadic grid storage
        c[j]=new double[N];                                             // scaling coefficients
        mask[j]=new bool[N];                                            // mask containing true for kept coefficients
    }									                                //
    for (j=0;j<J;j++) {                                                 //
        int N=jPnts(j)-1;                                               //
        d[j]=new double[N];                                             // detail coefficient
    }                                                                   //
    //------- Define timestep size -------------------------------------//
    int num_steps=40;                                            	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/num_steps;                                 	    // timestep size
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=jPnts(j-1)-1;                                             //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=2.*pow(2.,-(j+shift))*static_cast<double>(k);     // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Sample initial function on grid Gt -----------------------//
    for (k=0;k<jPnts(J);k++) c[J][k]=IC.f(x[J][k]);                	    // evaluate initial condition at collocation points 
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(x,c,d,J,interpPnts);                                      //
    //------- Remove coefficients below the threshold ------------------//
    for (j=0;j<J;j++) {                                                 //
        int N=jPnts(j);                                                 //
        for (k=0;k<N-1;k++) {                                           //
            if (abs(d[j][k])<threshold) {                       
                mask[j+1][2*k+1]=false;         // knock out points below threshold
                d[j][k]=0.;
            }
            else mask[j+1][2*k+1]=true;                                 // keep points above threshold, include in mask
        }                                                               //
    }                                                                   //
/*    //------- Include coarsest scaling coeff's in mask -----------------//
    for (j=0;j<=J;j++) {                                                 //
        int gridMultplr=pow(2,j-0);                                     //
        for (k=0;k<jPnts(0);k++) mask[j][gridMultplr*k]=true;           // all scaling coefficients at coarsest level included
    }                                                                   // */
/*    //------- Extend mask recursively ----------------------------------//
    for (j=J-2;j>0;j--) {
        int N=jPnts(j);
        for (k=0;k<N-1;k++) {
            if (mask[j+1][2*k+1]==true) {
                int leftPnt=-interpPnts+1+k;
                int rightPnt=interpPnts+k;
                while ( leftPnt < 0 ) {
                    leftPnt++;
                    rightPnt++;
                }
                while ( rightPnt > (N-1) ) {
                    leftPnt--;
                    rightPnt--;
                }
                for (int l=leftPnt;l<=rightPnt;l++) {
                    mask[j][k+l]=true;
                }
            }
        }
    } */
    //------- Calculate spatial derivatives ----------------------------//
    for (j=0;j<J;j++) {                                                 // 
        int N=jPnts(j);                                                 // number of points at current level
        int gridMultplr=pow(2,J-j);                                     // constant needed to get to same point at higher level
        for (k=0;k<N;k++) {                                             //
            if (mask[j][k]==true) {                                     // check if point in mask
                double xEval=x[j][k];                                   // evaluation point
                ux[gridMultplr*k]=lagrInterpD1(xEval,x[j],c[j],k,       // compute derivative from lagrange polynomial
                                    interpPnts,N);                      //
                activPnt[gridMultplr*k]=true;                           // represent this point at solution time
            }                                                           //
        }                                                               //
    }                                                                   //
    //------- Reconstruct function using wavelets ----------------------//    
    double** phi=new double*[J+1];
    double** psi=new double*[J+1];
    for (int j=0;j<=J;j++) {
        int N=jPnts(j);
        phi[j]=new double[N];
        psi[j]=new double[N];
    }
    for (j=0;j<J;j++) {
        int N=jPnts(j);
        for (int l=0;l<N;l++) {
            if (j==0) scaling_subd(phi,x,j,l,J,interpPnts);
            detail_subd(psi,x,j,l,J,interpPnts);
            for (i=0;i<jPnts(J);i++) {
                if (j==0) u[i]+=c[j][l]*phi[J][i];
                if (l<N-1) u[i]+=d[j][l]*psi[J][i]; 
            }
            phi[j][l]=0.;
            psi[j][l]=0.;
        }
    }
    delete[] phi;
    delete[] psi;
    //------- Output solution to file ----------------------------------//
    ofstream output1;                            	
    char fn1[30];                               		 
    snprintf(fn1,sizeof fn1,"output/solution.dat"); 			
    output1.open(fn1);                            	 
    for (int i=0;i<jPnts(J);i++) {  
        if (activPnt[i]==true) {
            output1<<x[J][i]<<" "<<u[i]<<endl;
        }
    }
    output1.close(); 
    //------- Output derivative to file -------------------------------//
    ofstream output2;                            	
    char fn[30];                               		 
    snprintf(fn,sizeof fn,"output/derivative.dat"); 			
    output2.open(fn);                            	 
    for (int i=0;i<jPnts(J);i++) {  
        if (activPnt[i]==true) {output2<<x[J][i]<<" "<<ux[i]<<endl;}
    }
    output2.close(); 
    //------- Output coefficient plot ---------------------------------//
    for (j=0;j<=J;j++) {
        ofstream output;
        char fn[25];
        snprintf(fn,sizeof fn,"output/coeff%d.dat",j);
        output.open(fn);
        if (j==0) {
            for (int i=0;i<jPnts(j);i++) {
                output<<x[j][i]<<" "<<0<<endl;
            }
        } else {
            for (int i=0;i<jPnts(j);i++) {
                if(mask[j][i]==true) {
                    output<<x[j][i]<<" "<<j<<endl;
                } else { 
                    output<<x[j][i]<<" "<<-10<<endl;
                }
            }
        }
        output.close();
    } 
    delete[] x;
    delete[] c;
    delete[] d;
    delete[] u;
    delete[] ux; 
    delete[] activPnt;
    return 0; 
}
