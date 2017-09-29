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

int shift;

// Inline functions:
int inline jPnts(int j) {return pow(2,j+shift)+1;}

int main(void) {
    //------- General parameters ---------------------------------------//
    //int numPoints=32;                                                 // number of level j=Jmax collocation points
    shift=2;
    int J=4;
    int interpPnts=2;                                                   // half the number of points used for interpolation
    //int J=log2(numPoints);                                           	// maximum scale level
    double threshold=pow(10.,-4);                               	    // error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    //------- Class definitions ----------------------------------------//
    initial_condition IC;                                           	// declare 'initial condition' class variable
    //------- Declare arrays -------------------------------------------//
    int N=jPnts(J);                                                     //
    double** u=new double*[J+1];                	                    // solution variable current timestep
    double** ux=new double*[J+1];                                       // first derivative of the solution wrt x
    double** x=new double*[J+1];			            		        // dyadic grid storage
    double** c=new double*[J+1];                                        // scaling function coefficients
    double** d=new double*[J];                                          // detail function coefficients
    bool** mask=new bool*[J];                                           // mask referencing all active grid points 
    for (j=0;j<=J;j++) {						                        //
    	int N=jPnts(j); 			    	                            // 
        u[j]=new double[N];                                             // solution variable old
        ux[j]=new double[N];                                            // first derivative
	    x[j]=new double[N];  						                    // dyadic grid storage
        c[j]=new double[N];                                             // scaling coefficients
        mask[j]=new bool[N];                                            // mask containing true for kept coefficients
    }									                                //
    for (j=0;j<J;j++) {                                                 //
        int N=jPnts(j)-1;                                               //
        d[j]=new double[N];                                             // detail coefficient
    }                                                                   //
    //------- Define timestep size -------------------------------------//
    int num_steps=1000;                                         	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/num_steps;                                 	    // timestep size
    //------- Populate dyadic grid -------------------------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=jPnts(j-1)-1;                                             //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=2.*pow(2.,-(j+shift))*k;                 	        // values of x on dyadic grid
            mask[j][k+N]=true;                                          //
        }                                                       	    //
    }                                                           	    //
    //------- Sample initial function on grid Gt -----------------------//
    for (k=0;k<N;k++) c[J][k]=IC.f(x[J][k]);                      	    // evaluate initial condition at collocation points 
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(x,u[J],c,d,J,interpPnts);                                 //
    //------- Remove coefficients below the threshold ------------------//
    for (j=0;j<J;j++) {                                                 
        int N=jPnts(j);                                               
        for (k=0;k<N;k++) {                                              
            if (abs(d[j][k])<threshold) mask[j+1][2*k+1]=false;
            else mask[j+1][2*k+1]=true;
        }
    }
    //------- Calculate first spatial derivative -----------------------//
    for (j=J-1;j>=0;j--) {
        int N=jPnts(j);
        for (k=0;k<N-1;k++) {
            if (mask[j+1][2*k+1]==false && (mask[j][k+1]==true || mask[j][k-1]==true)) {
                cout<<"j="<<j<<", k="<<k<<endl;
                double xEval=x[j+1][2*k+1];
                ux[j+1][2*k+1]=lagrInterpD1(xEval,x[j],c[j],k,interpPnts,N);
                xEval=x[j][k];
                ux[j+1][2*k]=lagrInterpD1(xEval,x[j],c[j],k,interpPnts,N);
            }
        }
    } 
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
                if (j==0) u[J][i]+=c[j][l]*phi[J][i];
                if (l<N-1) u[J][i]+=d[j][l]*psi[J][i]; 
            }
            phi[j][l]=0.;
            psi[j][l]=0.;
        }
    }
    delete[] phi;
    delete[] psi;
    //------- Output solution to file ---------------------------------//
    ofstream output;                            	
    char fn[30];                               		 
    snprintf(fn,sizeof fn,"output/solution.dat"); 			
    output.open(fn);                            	 
    for (int i=0;i<jPnts(J);i++) {  
        output<<x[J][i]<<" "<<ux[J][i]<<endl;     
    }
    output.close(); 
    //------- Output coefficient plot ---------------------------------//
    for (j=0;j<J;j++) {
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
