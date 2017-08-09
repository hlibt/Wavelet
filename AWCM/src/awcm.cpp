#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "wavelet.h"
#include "initial_conditions.h"
#include "awcm.hpp"
#define PI 3.14159265
using namespace std;

    //------------------------------------------------------------------------------//
    //                                                                              //
    //     ADAPTIVE WAVELET COLLOCATION METHOD USING 1-ST GENERATION WAVELETS TO    //
    //     SOLVE THE 1D VISCOUS BURGERS EQUATION ON AN ADAPTIVE-DYADIC GRID         //
    //                                                                              //
    //     AUTHOR: BRANDON GUSTO                                                    //
    //     DATE: Aug 01, 2017                                                       //
    //                                                                              //
    //------------------------------------------------------------------------------//

void time_stamp(int time,double diff,double dt);

int main(void) {
    //------- DEFINE GENERAL PARAMETERS --------------------------------//
    double v=0.0001;                                            	// kinematic diffusivity constant  
    int j;                                                      	// counter variable for wavelet level
    int k;                                                      	// counter variable for spatial index
    double threshold=5*pow(10.,-3);                             	// error tolerance for wavelet coefficients
    initial_condition IC;                                       	// declare 'initial condition' class variable
    wavelet db4;                                                	// declare 'wavelet' class variable
    //------- DECLARE ARRAYS -------------------------------------------//
    double sum;								// summation variable
    double sum1;							// summation variable
    double** residual=new double*[J+1];                         	// the residual between approximation Uj(x) and Uj-1(x) - create rows
    double** wave_coeff=new double*[J+1]				//
    double** U_old=new double*[J+1];					//
    double** Ux=new double*[J+1];					// first derivative of U with respect to x
    double** Uxx=new double*[J+1];					// second derivative of U with respect to x
    double** U_new=new double*[J+1];					// solution after time integration
    double** x=new double*[J+1];					// dyadic points
    double* a=new double[J+1];                                  	// wavelet dilation constants
    double** b=new double*[J+1];                                	// wavelet translation constants
    for (j=0;j<=J;j++) {						//
	int N=pow(2,j+2);						//
        U_old[j]=new double[N+1];					//
	U_new[j]=new double[N+1];					//
	Ux[j]=new double[N+1];						//
	Uxx[j]=new double[N+1];						//
	wave_coeff[j]=new double[N+1];					//
	residual[j]=new double[N+1];					//
	x[j]=new double[N+1];						//
	b[j]=new double[N+1]; 		                 		//
    }									//
    //------- DEFINE TIMESTEP SIZE -------------------------------------//
    int num_steps=1000;                                         	// number of timesteps     
    double ti=0.;                                               	// initial simulation time  
    double tf=1.;                                               	// final simulation time     
    double dt=(t_f-t_i)/num_steps;                              	// timestep size
    //------- SET UP DYADIC GRID ---------------------------------------//
    int num_points=64;                                                 	// number of level j=0 collocation points
    int J=log2(num_points);                                            	// number of scales possible
    int num_ext_wave_left=1;                                    	// number of external wavelets on left end of grid domain
    int num_ext_wave_right=1;                                   	// number of externam wavelets on right end of grid domain
    int L=1;                                                    	// indicates largest scale in wavelet basis
    double U_left_bound=0.;                                       	// left boundary point of the domain
    double U_right_bound=1.;                                      	// right boundary point of the domain
    double b0=1.;                                               	// wavelet translation constant
    double a0=pow(2.,-L)*(U_right_bound-U_left_bound)/b0;           	// wavelet dilation constant    
    for (j=0;j<=J;j++) {                                        	//
        a[j]=pow(2.,-(j+1.))*a0;                                	// 
        for (k=0;k<=pow(2,j+2);k++) {                           	//
            x[j][k]=U_left_bound+a[j]*b0*k;                       	// values of x on dyadic grid
        }                                                       	//
    }                                                           	//
    //------- SAMPLE INITIAL FUNCTION ON Gt ----------------------------//
    for (j=0;j<=J;j++) {                                        	//
        for (k=0;k<=pow(2,j+2);k++) {                           	//
            U_old[j][k]=IC.f(x[j][k]);                          	// evaluate initial condition at collocation points on dyadic grid
        }                                                       	//
    }                                                           	//
    //------- COMPUTE RESIDUALS ----------------------------------------//
    for (j=0;j<=J;j++) {                                        	// begin solving for coeff's one level at a time
        if (j==0) {                                             	//
            for (k=0;k<=pow(2,j+2);k++) {                       	//
                residual[j][k]=U_old[j][k];                     	// populate j=0 'rhs' of matrix system
            }                                                   	//
        }                                                       	//
        else if (j>0) {                                         	//
	    for (i=0;i<=pow(2,j+2);i++) {				// loop through all collocation points at current level
		sum=0.;							// create summation variable for contributions of lower levels
		for (int l=0;l<=j-1;l++) {				// loop from l=0 to l=j-1 (contributions of lower levels)
                    for (k=0;k<=pow(2,l+2);k++) {               	// loop through all spatial indices at level l 
			db4.set_params(l,k);				// set parameters to calculate daughter wavelet value
		        sum+=wave_coeff[l][k]*db4.daughter(x[j][i]);    // note that the operation is performed on grid level j 
                    }                                           	//
		}							//
		residual[j][i]=U_old[j][i]-sum; 			// populate j>0 'rhs' of matrix system
	    }								//
	}								//	
    //------- POPULATE WAVELET MATRIX Ajj ------------------------------//
	double** A=new double*[pow(2,j+2)+1];				// wavelet matrix filled with daughter wavelet functions
	for (i=0;i<=(pow(2,j+2));i++) A[i]=new double[pow(2,j+2)+1];	//
	for (i=0;i<=pow(2,j+2);i++) {					//
	    for (k=0;k<=pow(2,j+2);k++) {				//
		if (j>0) {						// when j>0
		    if (abs(wave_coeff[j][k])>threshold) {		// keep wavelet only if corresponding coefficient is above threshold
			db4.set_params(j,k);				// set parameters to calculate daughter wavelet
		    	A[i][k]=db4.daught(x[j][i]); 			// populate matrix A with daughter wavelet values
		    }							//
		    else {						//
		    	A[i][k]=0.;					// otherwise knock it out
		    }							//
		}							//
		else {							//
		    db4.set_params(j,k);				// set parameters to form daughter wavelet
		    A[i][k]=db4.daught(x[j][i],j,k);			// matrix A00 for when j=0;
		}							//
	    }								//
	}								//
    //------- SOLVE FOR WAVELET COEFFICIENTS AT j LEVEL ----------------//	
	wave_coeff[j]=BiCGSTAB(A,residual,tol,pow(2,j+2)+1,mxi);	// use the BiCGSTAB method to solve for wavelet coefficients
    //------- ELIMINATE COEFFICIENTS BELOW THRESHOLD -------------------//
   	for (k=0;k<=pow(2,j+2);k++) {					//
	    if (abs(wave_coeff[j][k])<=threshold) {			// check absolute value of coefficient against prescribed threshold
		wave_coeff[j][k]=0;					// make coefficient zero (for now - ideally would not have to store it)
//		if (j<J) wave_coeff[j+1][2*k]=0.;	       		// knock out coefficients at level j+1 (to make next matrix solve easier)
	    }								// 
	}								// 
    }									// end of j=0 to j=J loop
    //------- IF Gt != Gt+1, EVALUATE U AT NEW POINTS ------------------//
    for (j=0;j<=J;j++) {						//
	for (i=0;i<=pow(2,j+2);i++) {		         		//
	    sum=0.;							//
	    for (int l=0;l<j;l++) {					//
	    	for (k=0;k<=pow(2,l+2);k++) {				//
		    if (abs(wave_coeff[l][k])>threshold) {		// if coefficient is above threshold, the corresponding U must be calculated
		        sum+=wave_coeff[l][k]*db4.daughter(x[j][i],l,k);// 
		    }							//
		}							//
	    }								//
	    U_old[j][i]=sum;						//
	}								//
    }									//
    //------- COMPUTE DERIVATIVES USING BASIS --------------------------//
 //   for (j=0;j<=J;j++) {						//
//	for (i=0;i<=pow(2,j+2);i++) {		         		//
//	    sum=0.;							//
//	    sum1=0.;							//
//	    for (int l=0;l<j;l++) {					//
//	    	for (k=0;k<=pow(2,l+2);k++) {				//
//		    if (abs(wave_coeff[l][k])>threshold) {		// 
//		        sum+=wave_coeff[l][k]*psi_x_jk(x[j][i],l,k);	//
//			sum1+=wave_coeff[l][k]*psi_xx_jk(x[j][i],l,k);  // 
//		    }							//
//		}							//
//	    }								//
//	    Ux[j][i]=sum;						//
//	    Uxx[j][i]=sum1;						//
//	}								//
 //   }									//
/*    //------- INTEGRATE SOLUTION IN TIME -------------------------------//
    for (j=0;j<=J;j++) {						//
	for (i=0;i<=pow(2,j+2);i++) {					//
	    U_new[j][i]=U_old[j][i]+dt*(v*Uxx[j][i]-			//
			 U_old[j][i]*Ux[j][i]);				//
	}								//
    }									//
    //------- OUTPUT SOLUTION TO FILE ----------------------------------//    
    ofstream output;                            			//
    char fn[20];                               				// 
    snprintf(fn,sizeof fn,"../output/%04d.dat",s); 			//
    output.open(fn);                            			// 
    output<<0.<<" "<<bcic.f1(t[s])<<endl;       // left boundary                        //
            output<<x[l]<<" "<<U_old[l]<<endl;      // output solution to file              //
            Uxx_old[l]=Uxx_new[l];                  // update for next iteration            //
            Ux_old[l]=Ux_new[l];                    //                                      //
            U_old[l]=U_new[l];                      //                                      //
        output<<1.<<" "<<bcic.f2(t[s])<<endl;       //                                      //
        output.close();                             // close file                           //
        time_stamp(s,v,dt);                         // print info to screen                 //
    } */                                              // end of temporal iteration            //
    return 0;                                       //                                      //
}   

void time_stamp(int time,double diff,double dt) {
    cout << " " << endl;
    cout << "------------------------------"<<endl;
    cout << " kinematic diffusivity: " << diff << endl;
    cout << " current step: " << time << endl;
    cout << " time-step size: " << dt << endl;
    cout << " simulation time: " << dt*time << endl;
    cout << "------------------------------"<<endl;
    cout << " " << endl;
}
