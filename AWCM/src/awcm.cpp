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
double* complete_pivot(double** A,int size);
bool chk_conv(double** A,double* y,double* b,double tolerance,int n);
double inner_product(double* a, double* b,int n);
double* ADOTX(double** A,double* x,int n);
double L2norm(double* x,int n);

int main(void) {
    //------- DEFINE GENERAL PARAMETERS --------------------------------//
    double v=0.0001;                                                	// kinematic diffusivity constant  
    int num_points=16;                                                 	// number of level j=0 collocation points
    int J=log2(num_points);                                            	// number of scales possible
    int num_ext_wave_left=1;                                    	    // number of external wavelets on left end of grid domain
    int num_ext_wave_right=1;                                   	    // number of externam wavelets on right end of grid domain
    int L=1;                                                    	    // indicates largest scale in wavelet basis
    double U_left_bound=0.;                                       	    // left boundary point of the domain
    double U_right_bound=1.;                                      	    // right boundary point of the domain
    int i;                                                              // counter variable for spatial index
    int j;                                                          	// j is the counter variable for wavelet level
    int k;                                                          	// k is the counter variable for spatial index
    double threshold=5*pow(10.,-3);                                 	// error tolerance for wavelet coefficients
    double tol=pow(10.,-8.);                                            // BiCGSTAB convergence tolerance
    int mxi=5000;                                                       // maximum number of BiCGSTAB iterations
    initial_condition IC;                                           	// declare 'initial condition' class variable
    wavelet db4;                                                    	// declare 'wavelet' class variable
    //------- DECLARE ARRAYS -------------------------------------------//
    double sum;							                            	// summation variable
    double sum1;						                        	    // summation variable
    double** residual=new double*[J+1];                         	    // the residual between approximation Uj(x) and Uj-1(x) - create rows
    double** wave_coeff=new double*[J+1];    			                //
    double** U_old=new double*[J+1];					                //
    double** U_new=new double*[J+1];				            	    // solution after time integration
    double** x=new double*[J+1];			            		        // dyadic points
    for (j=0;j<=J;j++) {						                        //
	int N=pow(2,j+2);						                            //
        U_old[j]=new double[N+1];					                    //
	    U_new[j]=new double[N+1];					                    //
	    wave_coeff[j]=new double[N+1];					                //
	    residual[j]=new double[N+1];					                //
	    x[j]=new double[N+1];						                    //
    }									                                //
    //------- DEFINE TIMESTEP SIZE -------------------------------------//
    int num_steps=1000;                                         	    // number of timesteps     
    double ti=0.;                                               	    // initial simulation time  
    double tf=1.;                                               	    // final simulation time     
    double dt=(tf-ti)/num_steps;                                 	    // timestep size
    //------- POPULATE DYADIC GRID -------------------------------------//
    double b0=1.;                                               	    // wavelet translation constant
    double a0=pow(2.,-L)*(U_right_bound-U_left_bound)/b0;           	// wavelet dilation constant    
    double aj;                                                          //
    double bjk;                                                         // 
    for (j=0;j<=J;j++) {                                        	    //
        aj=pow(2.,-(j+1.))*a0;                                     	    // 
        for (k=0;k<=pow(2,j+2);k++) {                           	    //
            x[j][k]=U_left_bound+aj*b0*k;                         	    // values of x on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- SAMPLE INITIAL FUNCTION ON GRID 'Gt' ---------------------//
    for (j=0;j<=J;j++) {                                        	    //
        for (k=0;k<=pow(2,j+2);k++) {                           	    //
            U_old[j][k]=IC.f(x[j][k]);                          	    // evaluate initial condition at collocation points on dyadic grid
        }                                                       	    //
    }                                                           	    //
    //------- COMPUTE RESIDUALS ----------------------------------------//
    for (j=0;j<=J;j++) {                                        	    // begin solving for coeff's one level at a time
        int N=pow(2,j+2);                                               //
        if (j==0) {                                             	    //
            for (k=0;k<=N;k++) {                                   	    //
                residual[j][k]=U_old[j][k];                     	    // populate j=0 'rhs' of matrix system
            }                                                   	    //
        }                                                       	    //
        else if (j>0) {                                         	    //
	        for (i=0;i<=N;i++) {	        			                // loop through all collocation points at current level
		        sum=0.;							                        // summation variable for contributions of lower levels
		        for (int l=0;l<=j-1;l++) {				                // loop from l=0 to l=j-1 (contributions of lower levels)
                    for (k=0;k<=pow(2,l+2);k++) {                    	// loop through all spatial indices at level l 
			            db4.set_params(l,k);				            // set parameters to calculate daughter wavelet value
		                sum+=wave_coeff[l][k]*db4.daughter(x[j][i]);    // note that the operation is performed on grid level j 
                    }                                           	    //
		        }							                            //
		        residual[j][i]=U_old[j][i]-sum; 	                    // populate j>0 'rhs' of matrix system
	        }								                            //
	    }								                                //	
    //------- POPULATE WAVELET MATRIX Ajj ------------------------------//
	    double** A=new double*[N+1];        				            // wavelet matrix filled with daughter wavelet functions
	    for (i=0;i<=N;i++) {					                        //
            A[i]=new double[N+2];	                                    //
	        for (k=0;k<=N;k++) {		        		                //
		        db4.set_params(j,k);                                    // set parameters to calculate wavelets
		    	A[i][k]=db4.daughter(x[j][i]); 			                // populate matrix A with daughter wavelet values
	        }								                            //
	    }								                                //
        for (i=0;i<=N;i++) {                                            // augment matrix for guassian elimination algorithm
            A[i][N+1]=residual[j][i];                                   //
        }                                                               //
    //------- SOLVE FOR WAVELET COEFFICIENTS AT j LEVEL ----------------//	
	    wave_coeff[j]=complete_pivot(A,N+1);                            // wave_coeff[j]=BiCGSTAB(A,residual[j],tol,pow(2,j+2)+1,mxi);
    //------- ELIMINATE COEFFICIENTS BELOW THRESHOLD -------------------//
   	    for (k=0;k<=N;k++) {		        			                //
	        if (abs(wave_coeff[j][k])<=threshold) {			            // check absolute value of coefficient against prescribed threshold
                cout << j << k << endl;
		        wave_coeff[j][k]=0;					                    // make coefficient zero (for now - ideally would not have to store it)
		        if (j<J) {                                              //
                    wave_coeff[j+1][2*k]=0.;	       		            // knock out coefficients at level j+1 (to make next matrix solve easier)
                }								                        // 
	        }							                                // 
        }								                                // 
    }                                                                   // end of j=0 to j=J loop

    for (j=0;j<=J;j++) {
        for (i=0;i<=pow(2,j+2);i++) {                                                 	    //
            sum=0.;
            for (int l=0;l<=j;l++) {
                for (k=0;k<=pow(2,l+2);k++) {                           	    //
                    db4.set_params(l,k);
                    sum+=wave_coeff[l][k]*db4.daughter(x[j][i]);
                }
            }
            U_new[j][i]=sum;
        }              
    }                                                          	    //
    //------- IF Gt != Gt+1, EVALUATE U AT NEW POINTS ------------------//
  /*  for (j=0;j<=J;j++) {						//
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
    }					*/				//
    //------- COMPUTE DERIVATIVES USING BASIS --------------------------//
    //------- INTEGRATE SOLUTION IN TIME -------------------------------// 
    /*
    for (j=0;j<=J;j++) {						//
	for (i=0;i<=pow(2,j+2);i++) {					//
	    U_new[j][i]=U_old[j][i]+dt*(v*Uxx[j][i]-			//
			 U_old[j][i]*Ux[j][i]);				//
	}								//
    }									//
 */ //------- OUTPUT SOLUTION TO FILE ----------------------------------//    
    for (j=0;j<=J;j++) {
    ofstream output;                            			//
    char fn[20];                               				// 
    snprintf(fn,sizeof fn,"func%i.dat",j); 			//
    output.open(fn);                            			// 
    for (int l=0;l<=pow(2,j+2);l++) {
            output<<x[j][l]<<" "<<U_new[j][l]<<endl;      // output solution to file              //
    }
         output.close();                             // close file                           //     //  time_stamp(s,v,dt);                         // print info to screen                 //
    }
    return 0;                                       //                                      //
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

double* complete_pivot(double** A,int size) {       //                                      //
    int i, j, p, k, q, tmp;                         //                                      //
    int* NROW=new int[size];                        // pointer                              //
    int* NCOL=new int[size];                        // pointer                              //
    double** m=new double*[size];                   //                                      //
    for (i=0;i<size;i++) m[i]=new double[size+1];   //                                      //
    double* x=new double[size];                     // solution                             //
    double sigma;
    double A_pq;
    double A_kj;
    bool swap_row;
    bool swap_col;
    
    for (i=0;i<size;i++) {                          //                                      //    
        NROW[i]=i;                                  // initialize row pointer               //
        NCOL[i]=i;                                  // initialize column pointer            //
    }
    for (i=0;i<size;i++) {
        p=i;
        q=i;
        swap_row=false;
        swap_col=false;
        A_pq=abs(A[NROW[p]][NCOL[q]]);
        for (k=i;k<size;k++) {
            for (j=k;j<size;j++) {
                A_kj=abs(A[NROW[k]][NCOL[j]]);
                if (A_kj>A_pq) {
                    if (k!=p) {
                        p=k;
                        swap_row=true;
                    }
                    if (j!=q) {
                        q=j;
                        swap_col=true;
                    }
                    A_pq=abs(A[NROW[p]][NCOL[q]]);
                }
            }
        }
        if (swap_row==true) {
            tmp=NROW[i];
            NROW[i]=NROW[p];
            NROW[p]=tmp;
        }
        if (swap_col==true) {
            tmp=NCOL[i];
            NCOL[i]=NCOL[q];
            NCOL[q]=tmp;
        }
        for (j=i+1;j<size;j++) {
            m[NROW[j]][NCOL[i]]=A[NROW[j]][NCOL[i]]/A[NROW[i]][NCOL[i]];
            for (tmp=0;tmp<size+1;tmp++) {
                A[NROW[j]][tmp]=A[NROW[j]][tmp]-m[NROW[j]][NCOL[i]]*A[NROW[i]][tmp];
            }
        }
    }
    x[NCOL[size-1]]=A[NROW[size-1]][size]/A[NROW[size-1]][NCOL[size-1]];
    for (i=size-2;i>=0;i--) {
        sigma=0.;
        for (j=i;j<size;j++) {
            sigma+=A[NROW[i]][NCOL[j]]*x[NCOL[j]];
        }
        x[NCOL[i]]=(A[NROW[i]][size]-sigma)/A[NROW[i]][NCOL[i]];
    }
    return x;
}

   /* for (j=0;j<=J;j++) {
        for (i=0;i<=pow(2,j+2);i++) {                                                 	    //
            sum=0.;
            for (int l=0;l<=j;l++) {
                for (k=0;k<=pow(2,l+2);k++) {                           	    //
                    db4.set_params(l,k);
                    sum+=wave_coeff[l][k]*db4.daughter(x[j][i]);
                }
            }
            U_new[j][i]=sum;
        }              
    }  */                                                         	    //
