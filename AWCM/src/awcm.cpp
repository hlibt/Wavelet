#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "matrix_solvers/BiCGSTAB.h"
#include "conditions/conditions.hpp"
#include "transform/transform.hpp"
#include "interpolation/interpolation.hpp"
#include "phys_driver/phys_driver.hpp"
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
    //------- Grid and tolerance parameters ----------------------------//
    shift=2;                                                            // increases number of points of level j=0
    int J=10;                                                           // number of scales in the system
    int interpPnts=2;                                                   // half the number of points used for interpolation
    double threshold=pow(10.,-4);                               	    // error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j usually indicates decomposition scale
    int k;                                                          	// k is another variable for spatial index
    //------- Define timestep size -------------------------------------//
    int num_steps=40;                                            	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/num_steps;                                 	    // timestep size
    //------- Class declarations ---------------------------------------//
    initial_condition initCondition;                                    //
    //------- Declare arrays -------------------------------------------//
    int N=jPnts(J);                                                     // number of points at each level j
    double* u=new double[N];                      	                    // solution variable current timestep
    double* ux=new double[N];                                           // first derivative of the solution wrt x
    double* uxx=new double[N];                                          // second spatial derivative
    bool* activPnt=new bool[N];                                         // denotes points active on grid for output
    double** x=new double*[J+1];			            		        // dyadic grid storage (explicit storage)
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
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=jPnts(j-1)-1;                                             //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=2.*pow(2.,-(j+shift))*static_cast<double>(k);     // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Initially set c's to initial condition -------------------//
    for (k=0;k<jPnts(J);k++) c[J][k]=initCondition.f(x[J][k]);          // 
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(x,c,d,J,interpPnts);                                      // decompose signal into c's and d's
    //------- Remove coefficients below the threshold ------------------//
    thresholding(d,mask,threshold,J);                                   // knock out small d's
    //------- Include coarsest scaling coeff's in mask -----------------//
    for (k=0;k<jPnts(0);k++) mask[0][k]=true;                           // scaling coefficients at coarsest level included
    //------- Calculate spatial derivatives ----------------------------//
    for (j=0;j<J;j++) {                                                 // 
        int N=jPnts(j);                                                 // number of points at current level
        int gridMultplr=pow(2,J-j);                                     // constant needed to get to same point at higher level
        for (k=0;k<N;k++) {                                             //
            if (mask[j][k]==true) {                                     // check if point in mask
                double xEval=x[j][k];                                   // evaluation point
                ux[gridMultplr*k]=lagrInterpD1(xEval,x[j],c[j],k,       // compute derivative from lagrange polynomial
                                    interpPnts,N);                      //
                uxx[gridMultplr*k]=lagrInterpD2(xEval,x[j],c[j],k,       // compute derivative from lagrange polynomial
                                    interpPnts,N);                      //
                activPnt[gridMultplr*k]=true;                           // represent this point at solution time
            }                                                           //
        }                                                               //
    }                                                                   //
    //------- Reconstruct function using wavelets ----------------------//    
    reconstruction(x,u,c,d,J,interpPnts);                               // build the solution using wavelets
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
        if (activPnt[i]==true) {output2<<x[J][i]<<" "<<uxx[i]<<endl;}
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
