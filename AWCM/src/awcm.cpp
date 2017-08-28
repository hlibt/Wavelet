#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include "matrix_solvers/BiCGSTAB.h"
#include "initial_conditions/initial_conditions.hpp"
#include "wavelet_generation/wavelet_generation.hpp"
#include "transform/transform.hpp"
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
    int num_points=256;                                                 // number of level j=Jmax collocation points
    int J=log2(num_points);                                            	// maximum scale level
    double threshold=pow(10.,-2);                                 	// error tolerance for wavelet coefficients
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
    double** u_old=new double*[J+1];	            	                // solution at current/previous timestep
    double** u_new=new double*[J+1];				            	    // solution after time integration
    double** x=new double*[J+1];			            		        // dyadic points
    double** c=new double*[J+1];                                       // scaling function coefficients
    double** d=new double*[J+1];                                       // detail function coefficients
    for (j=0;j<=J;j++) {						                        //
	int N=pow(2,j+1);						                            // ** need to change this ***
        u_old[j]=new double[N+1];                                       //
        u_new[j]=new double[N+1];                                       //
	    x[j]=new double[N+1];						                    //
        c[j]=new double[N+1];                                           // scaling coefficients
        d[j]=new double[N+1];                                           // detail coefficients
    }									                                //
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j);                                                 //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=pow(2.,-j)*k;                             	    // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Sample initial function on grid Gt -----------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j);                                                 //
        for (k=-N;k<=N;k++) {                                   	    //
            u_old[j][k+N]=IC.f(x[j][k+N]);                        	    // evaluate initial condition at collocation points 
            u_new[j][k+N]=0.;
        }                                                       	    //
    }                                                           	    //
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(x,u_old[J],c,d,J,1);
    //------- Reconstruct function using wavelets ----------------------//    
    double** phi=new double*[J+1];
    double** psi=new double*[J+1];
    for (int i=0;i<=J;i++) {
        int N=pow(2,i+1)+1;
        phi[i]=new double[N];
        psi[i]=new double[N];
    }
    for (j=0;j<J;j++) {
        int N=pow(2,j+1)+1;
        for (int l=0;l<N;l++) {
            scaling_subd(phi,x,j,l,J,1);
            detail_subd(psi,x,j,l,J,1);
            for (i=0;i<=2*num_points;i++) {
                if (j==0) {
                    u_new[J][i]+=c[j][l]*phi[J][i]+d[j][l]*psi[J][i];
                } else {
                    u_new[J][i]+=d[j][l]*psi[J][i];
                }
            }
            phi[j][l]=0.;
            psi[j][l]=0.;
        }
    }
    delete[] phi;
    delete[] psi;
    //------- Output data to file --------------------------------------//
    ofstream output;                            	
    char fn[25];                               		 
    snprintf(fn,sizeof fn,"solution.dat"); 			
    output.open(fn);                            	 
    for (int t=0;t<=2*num_points;t++) {  
        output<<x[J][t]<<" "<<u_new[J][t]<<endl;     
    }
    output.close(); 
    delete[] x;
    delete[] c;
    delete[] d;
    delete[] u_old;
    delete[] u_new;
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
