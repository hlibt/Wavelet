#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include "haar.h"
#include "bc_ic.h"
#define PI 3.14159265
using namespace std;

    //----------------------------------------------//
    //                                              //
    //     1D VISCOUS BURGERS EQUATION SOLVER       //
    //     USING HAAR WAVELET COLLOCATION           //
    //                                              //
    //     AUTHOR: BRANDON GUSTO                    //
    //     DATE: JULY 13, 2017                      //
    //                                              //
    //----------------------------------------------//

double* BiCGSTAB(double** A,double* b,double tol,int size,int mxi);
double* ADOTX(double** A,double* b,int n);
double inner_product(double* a,double* b,int n);
bool chk_conv(double** A,double* y,double* b,double tolerance,int n);
double L2norm(double* x,int n);

int main(void) {
    //------- Declare initial parameters -----------//--------------------------------------//
    double v=0.1;                                   // kinematic diffusivity constant       //
    int M=16;                                        // half of collocation points           //
    int mxi=10000;                                    // max no. of iterations 4 matrix solvr //
    double tol=1.e-8;                               // matrix solver tolerance              //
    int J=log2(M);                                  // number of total scales               //
    int i;                                          // counter variable                     //
    int m;                                          //                                      //
    int l;                                          // counter variable                     //
    double c1,c2,c3,c4,c5;                          // constants                            //
    Haar H;                                         // declare class variable               //
    //------- Temporal discretization --------------//--------------------------------------//
    int N=100;                                       // number of timesteps                  //
    double t_i=0.;                                  // initial simulation time              //
    double t_f=1.;                                  // final simulation time                //   
    double dt=(t_f-t_i)/N;                          // timestep size                        //
    double* t=new double[N+1];                      // define temporal array                //
    for (i=0;i<=N;i++) t[i]=dt*i;                   // populate temporal array              //
    //------- Spatial discretization ---------------//--------------------------------------//
    double* x=new double[2*M];                      // array of collocation points          //
    for (i=1;i<=2*M;i++) x[i-1]=(i-.5)/(2.*M);      // populate array of collocation points //
    //------- Initialize matrices ------------------//--------------------------------------//
    double* c=new double[2*M];                      // array of wavelet coefficients        //
    double* b=new double[2*M];                      // RHS of matrix system                 //
    double** A=new double*[2*M];                    // initialize the LHS matrix 2Mx2M      //
    for (i=0;i<2*M;i++) A[i]=new double[2*M];       //                                      //
    double* U_old=new double [2*M];                 // initialize solution matrix U         //
    double* Ux_old=new double [2*M];                // initialize derivative of soln (U')   //
    double* Uxx_old=new double [2*M];               // initialize second derivtive of soln  //
    double* U_new=new double [2*M];                 // initialize solution matrix U         //
    double* Ux_new=new double [2*M];                // initialize derivative of soln (U')   //
    double* Uxx_new=new double [2*M];               // initialize second derivtive of soln  //
    //------- Initial conditions -------------------//--------------------------------------//
    BCIC bcic;                                      // declare class variable               //
    for (l=0;l<2*M;l++) {                           // populate the initial U,Ux,Uxx arrays //
        U_old[l]=bcic.f(x[l]);                      // initial function                     //
        Ux_old[l]=bcic.fx(x[l]);                    // derivative of intial function        //
        Uxx_old[l]=bcic.fxx(x[l]);                  // second derivative                    //
    }                                               //                                      //
    //------- Advance simulation in time -----------//--------------------------------------//
    for (int s=0;s<N;s++) {                         // march in time                        //
        for (l=0;l<2*M;l++) {                       // loop through collocation points      //
            c1=H.q1(x[l])-x[l]*H.q_tilda(1);        // and populate with functions for i=1  //
            c2=dt*(-v*H.h1(x[l])+Ux_old[l]*         //                                      //
                (H.q1(x[l])-x[l]*H.q_tilda(1)));    //                                      //
            c3=dt*U_old[l]*(H.p1(x[l])-             //                                      //
                H.q_tilda(1));                      //                                      //
            A[l][0]=c1+c2+c3;                       // populate LHS matrix for i=1          //
        }                                           //                                      //
        for (int j=0;j<=J;j++) {                    // repeat with wavelet functions (i!=1) //
            m=pow(2,j);                             //                                      //
            for (int k=0;k<m;k++) {                 // wavelet translation parameter        //
                H.set_params(k,m);                  // set parameters from k and m          //
                for (l=0;l<2*M;l++) {               // loop through collocation points      //
                    c1=H.q(x[l])-x[l]*H.q_tilda(m   //                                      //
                        +k+1);                      //                                      //
                    c2=dt*(-v*H.h(x[l])+Ux_old[l]*  //                                      //
                        (H.q(x[l])-x[l]*            //                                      //
                        H.q_tilda(m+k+1)));         //                                      //
                    c3=dt*U_old[l]*(H.p(x[l])-      //                                      //
                       H.q_tilda(m+k+1));           //                                      //
                    A[l][k+m]=c1+c2+c3;             // populate LHS matrix for i!=1         //
                }                                   //                                      //
            }                                       //                                      //
        }                                           //                                      //
        for (l=0;l<2*M;l++) {                       // populate RHS of matrix system        //
            c1=-bcic.f1t(t[s+1])-x[l]*              //                                      //
                (bcic.f2t(t[s+1])-bcic.f1t(t[s+1]));//                                      //
            c2=v*Uxx_old[l]-Ux_old[l]*(bcic.f1(t[s])//                                      //
                -bcic.f1(t[s+1]))+x[l]*(            //                                      //
                -bcic.f2(t[s+1])+bcic.f1(t[s+1])-   //                                      //
                 bcic.f2(t[s])+bcic.f1(t[s]))*      //                                      //
                 Ux_old[l];                         //                                      //
            c3=U_old[l]*(bcic.f1(t[s+1])+           //                                      //
                bcic.f1(t[s])-bcic.f2(t[s])-        //                                      //
                Ux_old[l]);                         //                                      //
            b[l]=c1+c2+c3;                          // RHS of matrix system                 //
        }                                           //                                      //
    //------- Solve matrix system ------------------//--------------------------------------//
        c=BiCGSTAB(A,b,tol,2*M,mxi);                // solve sys. for wavelet coefficients  //
    //------- Update solution variables ------------//--------------------------------------//
        for (l=0;l<2*M;l++) {                       // update solution variables            //
            c1=c[0]*H.h1(x[l]);                     // begin inner product for Uxx (scaling)//
            c2=c[0]*(H.p1(x[l])-H.q_tilda(1));      // begin inner product for Ux           //
            c3=c[0]*(H.q1(x[l])-x[l]*H.q_tilda(1)); // begin inner product for U            //
            for (int j=0;j<=J;j++) {                //                                      //
                m=pow(2,j);                         //                                      //
                for (int k=0;k<m;k++) {             //                                      //
                    i=m+k+1;                        //                                      //
                    H.set_params(k,m);              // set parameters Xi1, Xi2, Xi3         //
                    c1+=c[i-1]*H.h(x[l]);           // continue inner product calculation   //
                    c2+=c[i-1]*(H.p(x[l])-          // "                                    //
                         H.q_tilda(i));             // "                                    //
                    c3+=c[i-1]*(H.q(x[l])-x[l]*     // "                                    //
                         H.q_tilda(i));             // "                                    //
                }                                   //                                      //
            }                                       //                                      // 
            Uxx_new[l]=dt*c1+Uxx_old[l];            // update Uxx                           //
            Ux_new[l]=dt*c2+bcic.f2(t[s+1])-        //                                      //
                  bcic.f1(t[s+1])+bcic.f1(t[s])     //                                      //
                  -bcic.f2(t[s])+Ux_old[l];         //                                      //
            U_new[l]=dt*c3+x[l]*(bcic.f2(t[s+1])-   // update U                             //
                  bcic.f1(t[s+1])-bcic.f2(t[s])+    //                                      //
                  bcic.f1(t[s]))+bcic.f1(t[s+1])-   //                                      //
                  bcic.f1(t[s])+U_old[l];           //                                      //
        }                                           //                                      //
        ofstream output;                            //                                      //
        char fn[20];                                //                                      //
        snprintf(fn,sizeof fn,"output/%04d.dat",s); //                                      //
        output.open(fn);                            //                                      //
        for (l=0;l<2*M;l++) {                       //                                      //
            output<<x[l]<<" "<<U_old[l]<<endl;      //                                      //
            Uxx_old[l]=Uxx_new[l];                  // update for next iteration            //
            Ux_old[l]=Ux_new[l];                    //                                      //
            U_old[l]=U_new[l];                      //                                      //
        }                                           //                                      //
        output.close();                             // close file                           //
    }                                               // end of temporal iteration            //
    return 0; 
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
    double omega_new;
    double omega_old;
    double rho_new;
    double rho_old;
    double alpha;
    double beta;


    for (int i=0;i<size;i++) x0[i]=.0;              // populate initial solution guess
    double* c=ADOTX(A,x0,size);                     // calculate A*x0;
    for (int i=0;i<size;i++) r_old[i]=b[i]-c[i];    // start the residual with initial guess
    for (int i=0;i<size;i++) rhat[i]=r_old[i];      // choose arbitrary vector rhat
    rho_old=1.;
    omega_old=1.;
    alpha=1.;
    for (int i=0;i<size;i++) {
        v_old[i]=0.;
        p_old[i]=0.;
    }
    int k=1;
    bool iterate=true;
    do {
        rho_new=inner_product(rhat,r_old,size);
        beta=rho_new/rho_old*alpha/omega_old;
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

double L2norm(double* x,int n) {
    double sum=0.;
    for (int i=0;i<n;i++) sum+=sum+pow(x[i],2.);
    return sum=sqrt(sum);
}
