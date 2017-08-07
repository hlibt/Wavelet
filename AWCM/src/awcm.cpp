#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "wavelet.h"
#include "initial_conditions.h"
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
    //------- DEFINE GENERAL PARAMETERS ------------------------//
    double v=0.0001;                                            // kinematic diffusivity constant  
    int j;                                                      // counter variable for wavelet level
    int k;                                                      // counter variable for spatial index
    double threshold=5*pow(10.,-3);                             // error tolerance for wavelet coefficients
    initial_condition IC;                                       // declare 'initial condition' class variable
    wavelet db4;                                                // declare 'wavelet' class variable
    //------- DEFINE TIMESTEP SIZE -----------------------------//
    int num_steps=1000;                                         // number of timesteps     
    double ti=0.;                                               // initial simulation time  
    double tf=1.;                                               // final simulation time     
    double dt=(t_f-t_i)/num_steps;                              // timestep size
    //------- SET UP DYADIC GRID -------------------------------//
    int N=64;                                                   // number of level j=0 collocation points
    int J=log2(N);                                              // number of scales possible
    int num_ext_wave_left=1;                                    // number of external wavelets on left end of grid domain
    int num_ext_wave_right=1;                                   // number of externam wavelets on right end of grid domain
    int L=1;                                                    // indicates largest scale in wavelet basis
    double left_bound=0.;                                       // left boundary point of the domain
    double right_bound=1.;                                      // right boundary point of the domain
    double** x=new double*[J+1];                                // collocation points - create rows
    for (j=0;j<=J;j++) x[j]=new double[pow(2,j+2)+1];           // create columns
    double b0=1.;                                               // wavelet translation constant
    double a0=pow(2.,-L)*(right_bound-left_bound)/b0;           // wavelet dilation constant    
    double* a=new double[J+1];                                  // wavelet dilation constants
    double** b=new double*[J+1];                                // wavelet translation constants
    for (j=0;j<=J;j++) b[j]=new double[N+1];                    // fill rows of above ^^ with spatial points
    for (j=0;j<=J;j++) {                                        //
        a[j]=pow(2.,-(j+1.))*a0;                                // 
        for (k=0;k<=pow(2,j+2);k++) {                           //
            b[j][k]=left_bound+a[j]*b0*k;                       //  
            x[j][k]=b[j][k];                                    //
        }                                                       //
    }                                                           //
    //------- SAMPLE INITIAL FUNCTION ON Gt --------------------//
    double** U_old=new double*[J+1];                            // initialize solution matrix U - create rows
    for (j=0;j<=J;j++) U_old[j]=new double*[N];                 // create columns
    for (j=0;j<=J;j++) {                                        //
        for (k=0;k<=pow(2,j+2);k++) {                           //
            U_old[j][k]=IC.f(x[j][k]);                          // evaluate initial condition at collocation points on dyadic grid
        }                                                       //
    }                                                           //
    //------- COMPUTE THE RESIDUALS ----------------------------//
    double sum;							// summation variable
    double** residual=new double*[J+1];                         // the residual between approximation Uj(x) and Uj-1(x) - create rows
    for (j=0;j<=J;j++) residual[j]=new double[N];               // residual - create columns
    for (j=0;j<=J;j++) {                                        // begin solving for coeff's one level at a time
        if (j==0) {                                             //
            for (k=0;k<=pow(2,0+2);k++) {                       //
                residual[0][k]=U_old[0][k];                     // populate j=0 'rhs' of matrix system
            }                                                   //
        }                                                       //
        else if (j>0) {                                         //
	    for (i=0;i<=pow(2,j+2);i++) {			// loop through all collocation points at current level
		sum=0.;						//
                for (k=0;k<pow(2,j+1);k++) {                    //  
		    sum+=wave_coeff[j-1][k]			//
			*psi_jk(x[j][i],j-1,k);			//    	
                }                                               //
		residual[j][i]=residual[j-1][k]-sum;		// populate j>0 'rhs' of matrix system
	    }							//
	}							//	
    //------- POPULATE TRANSFORM MATRIX A ----------------------//
	double** A=new double*[pow(2,j+2)+1];			// wavelet transform matrix filled with daughter wavelet functions
	for (i=0;i<=(pow(2,j+2));i++) {				//
	    A[i]=new double[pow(2,j+2)+1];			//
	}							//
	for (i=0;i<=pow(2,j+2);i++) {				//
	    for (k=0;k<=pow(2,j+2);k++) {			//
		if (wave_coeff[j][k]>threshold) {		// keep wavelet only if corresponding coefficient is above threshold
		    A[i][k]=psi_jk(x[j][i],j,k);		//
		}						//
		else {						//
		    A[i][k]=0.;					// otherwise knock it out
		}						//
	    }							//
	}							//
    //------- SOLVE FOR WAVELET COEFFICIENTS AT j LEVEL --------//	
	double** wave_coeff=new double*[J+1]			//
	for (j=0;j<=J;j++) wave_coeff[j]=new double[pow(j+2)+1];//
	wave_coeff[j]=BiCGSTAB(A,residuals,tol,pow(2,j+2)+1,mxi);// use the BiCGSTAB method to solve for wavelet coefficients
    //------- ELIMINATE COEFFICIENTS BELOW THRESHOLD -----------//
   	for (k=0;k<=pow(2,j+2);k++) {				//
	    if (abs(wave_coeff[j][k])<=threshold) {		// check absolute value of coefficient against prescribed threshold
		wave_coeff[j+1][2*k]=0.;	        	// knock out coefficients at level j+1
	    }
	}
    //------- STEP : IF Gt != Gt+1, EVALUATE U AT NEW POINTS---//
    //------- STEP : INTEGRATE SOLUTION IN TIME ---------------//
    //double* U_new=new double[2*M];                              // initialize solution matrix U        
    //------- STEP : OUTPUT SOLUTION TO FILE ------------------//    
    ofstream output;                            //                                      //
    char fn[20];                                //                                      //
    snprintf(fn,sizeof fn,"../output/%04d.dat",s); //                                      //
    output.open(fn);                            //                                      //
    output<<0.<<" "<<bcic.f1(t[s])<<endl;       // left boundary                        //
            output<<x[l]<<" "<<U_old[l]<<endl;      // output solution to file              //
            Uxx_old[l]=Uxx_new[l];                  // update for next iteration            //
            Ux_old[l]=Ux_new[l];                    //                                      //
            U_old[l]=U_new[l];                      //                                      //
        output<<1.<<" "<<bcic.f2(t[s])<<endl;       //                                      //
        output.close();                             // close file                           //
        time_stamp(s,v,dt);                         // print info to screen                 //
    }                                               // end of temporal iteration            //
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
