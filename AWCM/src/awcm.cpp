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
    //------- STEP 1: DEFINE GENERAL PARAMETERS ----------------//
    double v=0.0001;                                            // kinematic diffusivity constant  
    int j;                                                      // counter variable for wavelet level
    int k;                                                      // counter variable for spatial index
    double eps=5*pow(10.,-3);                                   // error tolerance for wavelet coefficients
    initial_condition IC;                                       // declare 'initial condition' class variable
    wavelet db4;                                                // declare 'wavelet' class variable
    //------- STEP 2: DEFINE TIMESTEP SIZE ---------------------//
    int num_steps=1000;                                         // number of timesteps     
    double ti=0.;                                               // initial simulation time  
    double tf=1.;                                               // final simulation time     
    double dt=(t_f-t_i)/num_steps;                              // timestep size
    //------- STEP 3: SET UP DYADIC GRID -----------------------//
    int N=256;                                                  // number of level j=0 collocation points
    int J=log2(N);                                              // number of scales possible
    int num_ext_wave_left=1;                                    // number of external wavelets on left end of grid domain
    int num_ext_wave_right=1;                                   // number of externam wavelets on right end of grid domain
    int L=1;                                                    // indicates largest scale in wavelet basis
    double left_bound=0.;                                       // left boundary point of the domain
    double right_bound=1.;                                      // right boundary point of the domain
    double** x=new double*[J+1];                                // collocation points
    for (j=0;j<=J;j++) x[j]=new double[pow(2,j+2)+1];           //
    for (j=0;j<J;j++) {                                         //
        for (k=0;k<=pow(2,j+2);k++) {                           //
            x[j][k]=left_bound+k*pow(2,-j+2);                   // 
        }                                                       //
    }                                                           // 
    double b0=1.;                                               // wavelet translation constant
    double a0=pow(2.,-L)*(right_bound-left_bound)/b0;           // wavelet dilation constant    
    double* a=new double[J];                                    // wavelet dilation constants
    double** b=new double*[J+1];                                // wavelet translation constants
    for (j=0;j<=J;j++) b[j]=new double[N];                      // fill rows of above ^^ with spatial points
    for (j=0;j<=J;j++) {                                        //
        a[j]=pow(2.,-j)*a0;                                     // 
        for (k=0;k<pow(2,j);k++) {                              //
            b[j][k]=.5*(left_bound+right_bound)+a[j]*b0*k;      //  
        }                                                       //
    }                                                           //
    //------- STEP 4: SAMPLE INITIAL FUNCTION ON Gt ------------//
    double** U_old=new double*[J+1];                            // initialize solution matrix U        
    for (j=0;j<=J;j++) U_old[j]=new double*[N];                 //  
    for (j=0;j<=J;j++) {                                        //
        for (k=0;k<pow(2,j);k++) {                              //
            U_old[j][k]=IC.f(x[j][k]);                          //
        }                                                       //
    }                                                           //
    //------- STEP 5: COMPUTE THE RESIDUAL BETWEEN RESOLUTIONS -//
    double** residual=new double*[J+1];                         // the residual between approximation Uj(x) and Uj-1(x)
    for (j=0;j<=J;j++) residual[j]=new double*[N];              //
    for (j=0;j<=J;j++) {                                        //
        for (int l=0;l<=j;l++) {                                //
            for (k=0;k<pow(2.,j);k++) {                         //  
                if (j>=1 && l==j) {                             //
                    residual[j][k]=U_old[J][k]-U_old[j-1][k];   // equation 11 in Vasilyev & Paolucci 1997
                }                                               //
                else if (j>=1 && l<j) {                         // 
                    residual[j][k]=0.;                          // equation 11
                }                                               //
                else {                                          //
                    residual[j][k]=0.;                          //
                }                                               //
            }                                                   //
        }                                                       //
    }                                                           //
    double** residual=new double*[J];                           // the residual between approximation Uj(x) and Uj-1(x)
    for (j=0;j<J;j++) residual[j]=new double*[N];               //
    //------- STEP : ADJUST THE DYADIC GRID, Gt+1 -------------//    
    // adjust Gt+1 based on coefficients and epsilon parameter
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
