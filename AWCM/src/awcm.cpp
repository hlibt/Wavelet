#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "matrix_solvers/BiCGSTAB.h"
#include "initial_conditions/initial_conditions.hpp"
#include "wavelet_generation/scaling_subd.hpp"
#define PI 3.14159265
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

void time_stamp(int time,double diff,double dt);

int main(void) {
    //------- General parameters ---------------------------------------//
    int num_points=16;                                                 	// number of level j=0 collocation points
    int J=log2(num_points);                                            	// maximum scale level
    double u_bc1=-1.;                                            	    // left boundary point of the domain
    double u_bc2=1.;                                            	    // right boundary point of the domain
    double threshold=5*pow(10.,-3);                                 	// error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    initial_condition IC;                                           	// declare 'initial condition' class variable
    //------- Define timestep size -------------------------------------//
    int num_steps=1000;                                         	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/num_steps;                                 	    // timestep size
    //------- Declare arrays -------------------------------------------//
    double** scaling_coeff=new double*[J+1];    	                    // scaling wavelet coefficients
    double** detail_coeff=new double*[J];                               // detail wavelet coefficients
    double** u_old=new double*[J+1];	            	                // solution at current/previous timestep
    double** u_new=new double*[J+1];				            	    // solution after time integration
    double** x=new double*[J+1];			            		        // dyadic points
    int** map=new int*[J+1];                                            // pointer to deal with dynamic array indexing
    double* phi=new double[17];
    for (j=0;j<=J;j++) {						                        //
	int N=pow(2,j+1);						                            // ** need to change this ***
	    scaling_coeff[j]=new double[N+1];					            //
        if (j<J) detail_coeff[j]=new double[N+1];                       //
        u_old[j]=new double[N+1];                                       //
        u_new[j]=new double[N+1];                                       //
	    x[j]=new double[N+1];						                    //
    }									                                //
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j);                                                 //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=pow(2.,-j)*k;                                	    // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Sample initial function on grid Gt -----------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j);                                                 //
        for (k=-N;k<=N;k++) {                                   	    //
            u_old[j][k+N]=IC.f(x[j][k]);                          	    // evaluate initial condition at collocation points 
        }                                                       	    //
    }                                                           	    //
    phi=scaling_subd(x,1,1,4,2);
    ofstream output;                            	
    char fn[20];                               		 
    snprintf(fn,sizeof fn,"scaling.dat"); 			
    output.open(fn);                            	 
    for (int t=0;t<=num_points;t++) {
        output<<x[3][t]<<" "<<phi[t]<<endl;     
    }
    output.close();   
    return 0;
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
