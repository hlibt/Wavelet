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
#include "interpolation/interpolation.hpp"
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

int main(void) {
    //------- General parameters ---------------------------------------//
    int numPoints=256;                                                   // number of level j=Jmax collocation points
    int interpPnts=3;                                                   // half the number of points used for interpolation
    int J=log2(numPoints);                                            	// maximum scale level
    double threshold=pow(10.,-7);                               	    // error tolerance for wavelet coefficients
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    //------- Class definitions ----------------------------------------//
    initial_condition IC;                                           	// declare 'initial condition' class variable
    //------- Declare arrays -------------------------------------------//
    double** U=new double*[J+1];	            	                    // solution variable
    double** x=new double*[J+1];			            		        // dyadic grid storage
    double** c=new double*[J+1];                                        // scaling function coefficients
    double** d=new double*[J];                                          // detail function coefficients
    bool** mask=new bool*[J];                                           // mask referencing all active grid points 
    for (j=0;j<=J;j++) {						                        //
    	int N=pow(2,j+1)+1;			    	                            // 
        U[j]=new double[N];                                             // solution variable
	    x[j]=new double[N];  						                    // dyadic grid storage
        c[j]=new double[N];                                             // scaling coefficients
    }									                                //
    for (j=0;j<J;j++) {                                                 //
        int N=pow(2,j+1);                                               //
        d[j]=new double[N];                                             // detail coefficient
        mask[j]=new bool[N];                                            //
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
            x[j][k+N]=.5*pow(2.,-j)*k;                             	    // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Sample initial function on grid Gt -----------------------//
    for (j=0;j<=J;j++) {                                        	    //
        int N=pow(2,j+1)+1;                                             //
        for (k=0;k<N;k++) {                                        	    //
            U[j][k]=IC.f(x[j][k]);                                 	    // evaluate initial condition at collocation points 
        }                                                       	    //
    }                                                           	    //
    //------- Perform forward wavelet transform ------------------------//
    fwd_trans(x,U[J],c,d,J,interpPnts);                                 //
    //------- Remove coefficients below the threshold ------------------//
    for (j=0;j<J;j++) {                                                 //
        int N=pow(2,j+1);                                               //
        for (k=0;k<N;k++) {                                             // 
            if ( abs(d[j][k]) < threshold ) {
                mask[j][k]=false;
            } else {
                mask[j][k]=true;
            }
        }
    }
    //------- Calculate first spatial derivative -----------------------//
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
            if (j==0) {
                scaling_subd(phi,x,j,l,J,interpPnts);
            }
            detail_subd(psi,x,j,l,J,interpPnts);
            for (i=0;i<=2*numPoints;i++) {
                if (j==0) {
                    U[J][i]+=c[j][l]*phi[J][i];
                } 
                if (l<N-1) {
                    U[J][i]+=d[j][l]*psi[J][i];
                }
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
    snprintf(fn,sizeof fn,"solution.dat"); 			
    output.open(fn);                            	 
    for (int t=0;t<=2*numPoints;t++) {  
        output<<x[J][t]<<" "<<U[J][t]<<endl;     
    }
    output.close();
    //------- Output coefficient plot ---------------------------------//
    for (j=0;j<J;j++) {
        int N=pow(2,j+1);   
        ofstream output;
        char fn[25];
        snprintf(fn,sizeof fn,"coeff%d.dat",j);
        output.open(fn);
        if (j==0) {
            for (int t=0;t<N+1;t++) {
                output<<x[j][t]<<" "<<0<<endl;
            }
        } else {
            for (int t=0;t<N+1;t++) {
                if (t%2==1 && mask[j][t-1]==true) output<<x[j][t]<<" "<<j<<endl;
            }
        }
        output.close();
    }
    delete[] x;
    delete[] c;
    delete[] d;
    delete[] U;
    delete[] mask;
    return 0; 
}

/*    for (j=0;j<J;j++) {
        int n=pow(2,j+1)+1;
        for (k=0;k<n;k++) {
            if ( mask[j+1][2*k] == true ) {
            } else {
                double* xstar=new double[2*interpPnts];               
                double* fstar=new double[2*interpPnts];
                double searchWidth=0.0001;
                int numFound=0;
                do { 
                    numFound=0;
                    for (int jstar=J;jstar>=0;jstar--) {
                        int N=pow(2,jstar+1)+1;
                        for (int kstar=0;kstar<N;kstar++) {
                            double diff=abs(x[j][k]-x[jstar][kstar]);
                            if ( mask[jstar][kstar] == true && diff < searchWidth && numFound < 2*interpPnts ) {
                                xstar[numFound]=x[jstar][kstar];
                                fstar[numFound]=u_old[jstar][kstar];
                                numFound++;
                            }
                        }
                    }
                    searchWidth+=0.0001;
                } while ( numFound < 2*interpPnts ); // end while 
                cout<<"x is "<<x[j][k]<<endl;
                cout<<xstar[0]<<endl;
                cout<<xstar[1]<<endl;
                cout<<xstar[2]<<endl;
                cout<<xstar[3]<<endl;
                cout<<"end"<<endl; 
                // calculate derivative with lagrange
                delete[] xstar;
                delete[] fstar;
            }   
        }
    }  // end for  */
