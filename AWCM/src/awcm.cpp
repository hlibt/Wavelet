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

//function declarations
void time_stamp(int time,double diff,double dt);
double* BiCGSTAB(double** A,double* b,double tol,int size,int mxi);
bool chk_conv(double** A,double* y,double* b,double tolerance,int n);
double inner_product(double* a, double* b,int n);
double* ADOTX(double** A,double* x,int n);
double L2norm(double* x,int n);

int main(void) {
    //------- General parameters ---------------------------------------//
    int num_points=8;                                                 	// number of level j=0 collocation points
    int J=log2(num_points);                                            	// maximum scale level
    double u_bc1=-1.;                                            	    // left boundary point of the domain
    double u_bc2=1.;                                            	    // right boundary point of the domain
    double threshold=5*pow(10.,-3);                                 	// error tolerance for wavelet coefficients
    double even_weight=0.5;                                             // interpolating coefficients (lagrange poly.)
    double odd_weight=0.5;                                              //
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    initial_condition IC;                                           	// declare 'initial condition' class variable
    wavelet db4;                                                    	// declare 'wavelet' class variable
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
        N=pow(2,j);                                                     //
        for (k=-N;k<=N;k++) {                                    	    //
            x[j][k+N]=pow(2.,j)*k;                                	    // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- Sample the scaling function phi with transform -----------//
    for (j=0;j<=J;j++) {                                                //
        int N=pow(2,j);                                                 //
        for (int m=-N;m<=N;m++) {                                       //
            for (k=-N;k<=N;k++) {                                       //
                scaling_coeff[j][k+N]=kronicker_delta(k,m);             // set the scaling coefficients to the kronicker delta
            }                                                           //
            for (int jp=j;jp<=J-1;jp++) {                               // now perform the inverse transform to get the phi_jm's
                int Np=pow(2,jp);                                       //
                for (k=-Np;k<=Np;k++) {                                 //
                    scaling_coeff[jp+1][2*k+2*Np]=scaling_coeff[jp][k]; //
                    double tmp=0.;                                      //
                    for (int l=-neighbors+1;l<=neighbors;l++) {         //
                        tmp+=even_weight*scaling_coeff[j][k+l+Np];      //
                    }                                                   //
                    scaling_coeff[jp+1][2*k+2*Np+1]=tmp;                //
                }                                                       //
                phi[j][m]=scaling_coeff[j][m];                          //
            }                                                           //
        }                                                               //
    }
    //------- Sample initial function on grid Gt -----------------------//
    for (j=0;j<=J;j++) {                                        	    //
        N=pow(2,j);                                                     //
        for (k=-N;k<=N;k++) {                                   	    //
            u_old[j][k+N]=IC.f(x[j][k]);                          	    // evaluate initial condition at collocation points on dyadic grid
        }                                                       	    //
    }                                                           	    //
}

double* BiCGSTAB(double** A,double* b,double tol,int size,int mxi) {
    // Purpose:  Solves the matrix equation 'Ax=b' using
    //          the Biconjugate gradient stabilized method. 
    // Modified: July 03, 2017
    //
    // Author: Brandon Gusto
    //
    // Parameters: Input --> Diag, the elements on the diagonal of the matrix 'A'
    //             Input --> Lowr, the elements on the lower diagonal of the matrix 'A'
    //             Input --> Uppr, the elements on the upper diagonal of the matrix 'A'
    //             Input --> b, the elements in the vector 'b' in the RHS of 'Ax=b'
    //             Input --> x, the solution vector
    //             Input --> tol, the tolerance with which the BiCGSTAB method solves the system
    //             Input --> mxi, the maximum number of iterations
    //             Local -->
    double* r_new=new double[size];
    double* r_old=new double[size];
    double* rhat=new double[size];
    double* v_new=new double[size];    
    double* v_old=new double[size];
    double* p_new=new double[size];
    double* p_old=new double[size];
    double* h=new double[size];
    double* s=new double[size];
    double* t=new double[size];
    double* x0=new double[size];
    double* tmp=new double[size];
    double omega_new;
    double omega_old;
    double rho_new;
    double rho_old;
    double alpha;
    double beta;

    for (int i=0;i<size;i++) x0[i]=10.0;             // populate initial solution guess
    tmp=ADOTX(A,x0,size);                           // calculate A*x0;
    for (int i=0;i<size;i++) r_old[i]=b[i]-tmp[i];  // start the residual with initial guess
    for (int i=0;i<size;i++) rhat[i]=r_old[i];      // choose arbitrary vector rhat
    rho_old=1.; omega_old=1.; alpha=1.;             // initial parameters
    for (int i=0;i<size;i++) {
        v_old[i]=0.;
        p_old[i]=0.;
    }
    int k=1;
    bool iterate=true;
    do {
        rho_new=inner_product(rhat,r_old,size);
        beta=(rho_new/rho_old)*(alpha/omega_old);
        for (int i=0;i<size;i++) p_new[i]=r_old[i]+beta*(p_old[i]-omega_old*v_old[i]);
        v_new=ADOTX(A,p_new,size);
        alpha=rho_new/inner_product(rhat,v_new,size);
        for (int i=0;i<size;i++) h[i]=x0[i]+alpha*p_new[i];
        if (chk_conv(A,h,b,tol,size)==true) {
            iterate=false;
            for (int i=0;i<size;i++) x0[i]=h[i];
            cout << "Check 1 \n";
        }
        else {
            for (int i=0;i<size;i++) s[i]=r_old[i]-alpha*v_new[i];
            t=ADOTX(A,s,size);
            omega_new=inner_product(t,s,size)/inner_product(t,t,size);
            for (int i=0;i<size;i++) x0[i]=h[i]+omega_new*s[i];
            if (chk_conv(A,x0,b,tol,size)==true) {
                iterate=false;
                cout << "Check 2 \n";
            }
            else {
                for (int i=0;i<size;i++) r_new[i]=s[i]-omega_new*t[i];
            }
        }
        for (int i=0;i<size;i++) {
            r_old[i]=r_new[i];
            p_old[i]=p_new[i];
            v_old[i]=v_new[i];
        }            
        rho_old=rho_new;
        omega_old=omega_new;
        k++;
        if (k>mxi) {
            iterate=false;
            cout << "Max iterations reached \n";
        }
    }while(iterate==true);
    return x0;
}

double* ADOTX(double** A,double* x,int n) {
    double* output=new double[n];
    double sum=0.;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            sum+=A[i][j]*x[j];
        }
        output[i]=sum;
        sum=0;
    }
    return output;
}
 
double inner_product(double* a, double* b,int n) {
    double output=0.;
    for (int i=0;i<n;i++) output+=a[i]*b[i];
    return output;
}

double L2norm(double* x,int n) {
    double sum=0.;
    for (int i=0;i<n;i++) sum+=sum+pow(x[i],2.);
    return sum=sqrt(sum);
    for (int i=0;i<n;i++) sum+=pow(x[i],2.);
    return sqrt(sum);
}

bool chk_conv(double** A,double* y,double* b,double tolerance,int n) {
    double* residual=new double[n];
    double rel_error;
    bool conv;
    residual=ADOTX(A,y,n);
    for (int i=0;i<n;i++) residual[i]=b[i]-residual[i];
    rel_error=L2norm(residual,n)/L2norm(b,n);
    if (rel_error<tolerance)
        conv=true;
    else
        conv=false;
    return conv;
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

double kronicker_delta(int k, int m) {
    if (k==m) {
        return 1.;
    } else {
        return 0.;
    }
}

double scaling_subd(int j,int m,int k,int Jmax,npnts) {
    
    //------------------------------------------------------------------//
    // Information: scaling_subd performs the interpolating subdivision algorithm
    //              in order to determine the scaling functions phi_j,m sampled at the 
    //              specific locations of x_Jmax,k. 
    //
    // Input: 
    //              j     - level of the scaling function
    //              m     - translaton parameter of the scaling function
    //              Jmax  - maximum desired grid level for the point x
    //              k     - spatial index of x
    //              npnts - half the number of nearest points to use in subdivision scheme
    // Output:
    //              phi_j,m(x_Jmax,k)
    //------------------------------------------------------------------//     
                                                                        //
    double** c=new double*[Jmax];                                       // coefficients for interpolating subdivision
    for (int i=j;i<=Jmax;i++) {                                         //
        int n=pow(2,i+1)+1;                                             // number of points at level j
        c[i]=new double[n];                                             // intialize columns of c
    }                                                                   //
    int n=pow(2,j+1)+1;                                                 //
    for (int i=0;i<n;i++) {                                             //
        c[j][i]=kronicker_delta(k,m);                                   // set coefficients to kronicker delta function
    }                                                                   //
    for (int jstar=j;jstar<=Jmax;jstar++) {                             // begin inverse transform process
        int n=pow(2,jstar+1)+1;                                         // number of points at level jstar
        for (int i=0;i<n;i++) {                                         // 
            c[jstar+1][2*i]=c[jstar][i];                                // even points stay the same
            double tmp=0.;                                              // summation variable
            for (int l=-npnts+1;l<=npnts;l++) {                         // 
                tmp+=weight*c[jstar][i+l];                              //
            }                                                           //
            c[jstar+1][2*i+1]=tmp;                                      // odd points
        }                                                               //
    }                                                                   //
    return c[Jmax];                                                     // the final scaling function at sampled points
}                                                                       //
