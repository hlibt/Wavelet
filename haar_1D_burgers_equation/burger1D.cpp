#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <cmath>
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

double* complete_pivot(double** A,int size);

int main(void) {
    //------- Declare initial parameters -----------//--------------------------------------//
    double v=0.01;                                   // kinematic diffusivity constant       //
    int M=128;                                        // half of collocation points           //
    int J=log2(M);                                  // number of total scales               //
    int i;                                          // counter variable                     //
    int m;                                          //                                      //
    int l;                                          // counter variable                     //
    double c1,c2,c3;                                // constants                            //
    Haar H;                                         // declare class variable               //
    //------- Temporal discretization --------------//--------------------------------------//
    int N=1000;                                       // number of timesteps                  //
    double t_i=0.;                                  // initial simulation time              //
    double t_f=3.;                                  // final simulation time                //   
    double dt=(t_f-t_i)/N;                          // timestep size                        //
    double* t=new double[N+1];                      // define temporal array                //
    for (i=0;i<=N;i++) t[i]=dt*i;                   // populate temporal array              //
    //------- Spatial discretization ---------------//--------------------------------------//
    double* x=new double[2*M];                      // array of collocation points          //
    for (i=1;i<=2*M;i++) x[i-1]=(i-.5)/(2.*M);      // populate array of collocation points //
    //------- Initialize matrices ------------------//--------------------------------------//
    double* c=new double[2*M];                      // array of wavelet coefficients        //
    double** A=new double*[2*M];                    // initialize the LHS matrix 2Mx2M      //
    for (i=0;i<2*M;i++) A[i]=new double[2*M+1];     //                                      //
    double* U_old=new double[2*M];                  // initialize solution matrix U         //
    double* Ux_old=new double[2*M];                 // initialize derivative of soln (U')   //
    double* Uxx_old=new double[2*M];                // initialize second derivtive of soln  //
    double* U_new=new double[2*M];                  // initialize solution matrix U         //
    double* Ux_new=new double[2*M];                 // initialize derivative of soln (U')   //
    double* Uxx_new=new double[2*M];                // initialize second derivtive of soln  //
    //------- Initial conditions -------------------//--------------------------------------//
    BCIC bcic;                                      // declare class variable               //
    for (l=0;l<2*M;l++) {                           // populate the initial U,Ux,Uxx arrays //
        U_old[l]=bcic.f(x[l]);                      // initial function                     //
        Ux_old[l]=bcic.fx(x[l]);                    // derivative of intial function        //
        Uxx_old[l]=bcic.fxx(x[l]);                  // second derivative                    //
    }                                               //                                      //
    //------- Advance simulation in time -----------//--------------------------------------//
    for (int s=0;s<N;s++) {                         // march in time                        //
   cout<<s<<endl;
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
            A[l][2*M]=c1+c2+c3;                     // RHS of matrix system                 //
        }                                           //                                      //
    //------- Solve matrix system ------------------//--------------------------------------//
        c=complete_pivot(A,2*M);                // solve sys. for wavelet coefficients  //
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

double* complete_pivot(double** A,int size) {
    int i, j, p, k, q, tmp;
    int* NROW=new int[size];
    int* NCOL=new int[size];
    double** m=new double*[size];
    for (i=0;i<size;i++) m[i]=new double[size+1];
    double* x=new double[size];
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
