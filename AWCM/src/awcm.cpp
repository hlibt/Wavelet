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
void diffWave(double** f,double** df,double h,int J,int N);

    //------------------------------------------------------------------------------//
    //                                                                              //
    //     ADAPTIVE WAVELET COLLOCATION METHOD USING 2ND GENERATION WAVELETS TO     //
    //     SOLVE THE 1D VISCOUS BURGERS EQUATION ON AN ADAPTIVE-DYADIC GRID         //
    //                                                                              //
    //     AUTHOR: BRANDON GUSTO                                                    //
    //     DATE: Aug 01, 2017                                                       //
    //                                                                              //
    //------------------------------------------------------------------------------//

int main(void) {
    //------- General parameters ---------------------------------------//
    int numPoints=2;                                                    // number of level j=Jmax collocation points
    int interpPnts=1;                                                   // half the number of points used for interpolation
    int J=log2(numPoints);                                           	// maximum scale level
    double threshold=pow(10.,-4);                               	    // error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    //------- Class definitions ----------------------------------------//
    initial_condition IC;                                           	// declare 'initial condition' class variable
    //------- Declare arrays -------------------------------------------//
    double* Du1=new double[2*numPoints];                                // first derivative of the solution
    double* u_new=new double[2*numPoints];                              // solution at t+dt timestep
    double** u_old=new double*[J+1];            	                    // solution variable current timestep
    double** x=new double*[J+1];			            		        // dyadic grid storage
    double** c=new double*[J+1];                                        // scaling function coefficients
    double** d=new double*[J];                                          // detail function coefficients
    bool** mask=new bool*[J];                                           // mask referencing all active grid points 
    for (j=0;j<=J;j++) {						                        //
    	int N=pow(2,j+1)+1;			    	                            // 
        u_old[j]=new double[N];                                         // solution variable old
	    x[j]=new double[N];  						                    // dyadic grid storage
        c[j]=new double[N];                                             // scaling coefficients
        mask[j]=new bool[N];                                            // mask containing true for kept coefficients
    }									                                //
    for (j=0;j<J;j++) {                                                 //
        int N=pow(2,j+1);                                               //
        d[j]=new double[N];                                             // detail coefficient
    }                                                                   //
    //------- Define timestep size -------------------------------------//
    int num_steps=1000;                                         	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/num_steps;                                 	    // timestep size
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j);                                                 //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=pow(2.,-j)*k;                             	    // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Sample initial function on grid Gt -----------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j+1)+1;                                             //
        for (k=0;k<N;k++) {                                        	    //
            u_old[j][k]=IC.f(x[j][k]);                             	    // evaluate initial condition at collocation points 
            mask[j][k]=true;                                            // use loop also to initialize mask
        }                                                       	    //
    }                                                           	    //
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(x,u_old[J],c,d,J,interpPnts);                             //
    //------- Remove coefficients below the threshold ------------------//
    for (j=0;j<J;j++) {                                                 
        int N=pow(2,j+1);                                               
        for (k=0;k<N;k++) {                                              
            if ( abs(d[j][k]) < threshold ) {
                mask[j+1][2*k+1]=false;
            } else {
                mask[j+1][2*k+1]=true;
            }
        }
    }
    //------- Calculate first spatial derivative -----------------------//
    for (j=J-1;j>=0;j--) {
        int N=pow(2,j+1)+1;
        for (k=0;k<N-1;k++) {
            if (mask[j+1][2*k+1]==false) {
                double xEval=x[j+1][2*k+1];
                Du1[k]=lagrInterpD1(xEval,x[j],c[j],k,interpPnts,N);
		        cout<<"x is "<<x[j+1][2*k+1]<<" Du is "<<Du1[k]<<endl;
            }
        }
    }
    //------- Reconstruct function using wavelets ----------------------//    
    double** phi=new double*[J+1];
    double** psi=new double*[J+1];
    double** Dphi=new double*[J+1];
    double** Dpsi=new double*[J+1];
    for (int i=0;i<=J;i++) {
        int N=pow(2,i+1)+1;
        phi[i]=new double[N];
        psi[i]=new double[N];
        Dphi[i]=new double[N];
        Dpsi[i]=new double[N];
    }
    for (j=0;j<J;j++) {
        int N=pow(2,j+1)+1;
        for (int l=0;l<N;l++) {
            if (j==0) {
                scaling_subd(phi,x,j,l,J,interpPnts);
            }
            detail_subd(psi,x,j,l,J,interpPnts);
            diffWave(phi,Dphi,abs(x[J][1]-x[J][0]),J,pow(2,J+1)+1);
            diffWave(psi,Dpsi,abs(x[J][1]-x[J][0]),J,pow(2,J+1)+1);
            for (i=0;i<=2*numPoints;i++) {
                if (j==0) {
                    u_new[i]+=c[j][l]*phi[J][i];
                } 
                if (l<N-1) {
                    u_new[i]+=d[j][l]*psi[J][i];
                }
            }
            phi[j][l]=0.;
            psi[j][l]=0.;
            Dphi[j][l]=0.;
            Dpsi[j][l]=0.;
        }
    }
    delete[] phi;
    delete[] psi;
    delete[] Dphi;
    delete[] Dpsi;
    //------- Output solution to file ---------------------------------//
    ofstream output;                            	
    char fn[30];                               		 
    snprintf(fn,sizeof fn,"solution.dat"); 			
    output.open(fn);                            	 
    for (int t=0;t<=2*numPoints;t++) {  
        output<<x[J][t]<<" "<<u_new[t]<<endl;     
    }
    output.close();
    //------- Output coefficient plot ---------------------------------//
    for (j=0;j<J;j++) {
        int N=pow(2,j);   
        ofstream output;
        char fn[25];
        snprintf(fn,sizeof fn,"coeff%d.dat",j);
        output.open(fn);
        if (j==0) {
            for (int t=0;t<pow(2,j+1)+1;t++) {
                output<<x[j][t]<<" "<<0<<endl;
            }
        } else {
            for (int t=0;t<N;t++) {
                if(mask[j-1][t]==true) {
                    output<<x[j][2*t+1]<<" "<<j<<endl;
                } else { 
                    output<<x[j][2*t+1]<<" "<<-10<<endl;
                }
            }
        }
        output.close();
    }
    delete[] x;
    delete[] c;
    delete[] d;
    delete[] u_old;
    delete[] u_new;
    delete[] mask;
    delete[] Du1;
    return 0; 
}

void diffWave(double** f,double** df,double h,int J,int N) {
    df[J][0]=(-3.*f[J][0]+4.*f[J][1]-f[J][2])/(2.*h);
    for (int i=1;i<N-1;i++) {
        df[J][i]=(f[J][i+1]-f[J][i-1])/(2.*h);
    }
    df[J][N-1]=(3.*f[J][N-1]-4.*f[J][N-2]+f[J][N-2])/(2.*h);
    return;
}
