#include <cmath>
#include <string>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
#include "../interpolation/interpolation.hpp"

double rhs(double* gridpts, double* funcpts, int nactive, int i, double mod, double alpha, double c, string equation) { 

    //--------------------------------------------------------------------------------------//
    // Information: rhs.cpp computes the right-hand-side of the initial value problem
    //              at each timestep for a single grid point.                  
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.x               - the x-location of the collocation point
    //              collPnt.ux              - the derivative of 'u' to be computed 
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------------------//     

<<<<<<< HEAD
    //------- initialize vectors ---------------------------------------//
    double* gridPnts = new double[2*interpPnts];                        // the interpolating stencil
    double* funcPnts = new double[2*interpPnts];                        // the function evaluated at stencil points

    //------- determine if given stencil is adequate -------------------//
    if ( j!=J && collPnt[j+1][2*i+1].isMask == true 
            && ) {                 //
        j++;                                                            // increment wavelet-grid level
        i*=2;                                                           // represent grid index at next level
    }

    //------- modify the gridpoint's stencil values --------------------//
    double f0 = collPnt[j][i].u + mod;                                  // solution at time subinterval
    int p = (i-1)/2 - interpPnts + 1;                                   // starting grid index
    for (int k=0;k<2*interpPnts;k++) {                                  // loop through stencil
        gridPnts[k] = collPnt[j-1][p].x;                                //
        funcPnts[k] = collPnt[j-1][p].u + mod;                          // add modification to solution variable
=======
    //------- compute interpolation stencil ----------------------------//
    int N = 2*interpPnts + 1;                                           // number of points in the stencil
    double* xstencil = new double[N];                                   // the interpolating stencil
    double* fstencil = new double[N];                                   // the function value at stencil points
    double xeval = gridpts[i];                                          // point to evaluate derivative at
    double f0 = funcpts[i] + mod;                                       // solution at time subinterval
    int p = i - interpPnts;                                             // starting grid index
    while ( p < 0 ) { p++; }                                            // this loop checks to ensure the stencil consists of points in the domain
    while ( p + 2*interpPnts > (nactive-1) ) { p--; }                   // adjust stencil one point to the left if necessary
    for (int k=0;k<=2*interpPnts;k++) {                                 // populate stencil
        xstencil[k] = gridpts[p];                                       //
        fstencil[k] = funcpts[p] + mod;                                 //
>>>>>>> 2f6e487d56af862323ff001ca16a8666c02e1e39
        p++;                                                            //
    }                                                                   //

    //------- compute first-order spatial derivative -------------------//
    double f1 = 0.;                                                     // the first derivative of u
    for (p=0;p<N;p++) {                                                 // loop through stencil
        double sum = 0.;                                                //
        for (int k=0;k<N;k++) {                                         //
            if (k!=p) {                                                 //
                double productUpper = 1.;                               //
                for (int l=0;l<N;l++) {                                 //
                    if ( l!=k && l!=p ) {                               //
                        productUpper *= ( xeval - xstencil[l] );        //
                    }                                                   //
                }                                                       //
                sum += productUpper;                                    //
            }                                                           //
        }                                                               //
        double productLower = 1.;                                       //
        for (int k=0;k<N;k++) {                                         //
            if (k!=p) {                                                 //
                productLower *= ( xstencil[p] - xstencil[k] );          //
            }                                                           //
        }                                                               //
        f1 += fstencil[p] * sum / productLower;                         //
    }                                                                   //

    //------- compute second-order spatial derivative ------------------//
    double f2 = 0.;                                                     // the second derivative of u
    for (p=0;p<N;p++) {                                                 //
        double sumOuter = 0.;                                           //
        for (int l=0;l<N;l++) {                                         //
            if (l!=i) {                                                 //
                double sumInner = 0.;                                   //
                for (int m=0;m<N;m++) {                                 //
                    if (m!=i && m!=l) {                                 //
                        double product=1.;                              //
                        for (int k=0;k<N;k++) {                         //
                            if (k!=i && k!=l && k!=m) {                 //
                                product *= ( xeval - xstencil[k] ) /    //
                                          ( xstencil[i] - xstencil[k]); //
                            }                                           //
                        }                                               //
                        sumInner += product / (xstencil[i]-xstencil[m]);//
                    }                                                   //
                }                                                       //
                sumOuter += sumInner / ( xstencil[i]-xstencil[l] );     //
            }                                                           //
        }                                                               //
        f2 += fstencil[i] * sumOuter;                                   //
    }                                                                   //

    //------- Combine derivatives to form rhs --------------------------//
    if ( equation.compare(0,7,"burgers") == 0 ) {
        return -f0*f1 + alpha*f2;
    } else if ( equation.compare(0,16,"modified_burgers") == 0 ) {
        return ( -(f0 + c)*f1 + alpha*f2 );
    } else if ( equation.compare(0,19,"advection_diffusion") == 0 ) {
        return -c*f1 + alpha*f2;
    } else if ( equation.compare(0,9,"advection") == 0 ) {
        return -c*f1;
    } else if ( equation.compare(0,9,"diffusion") == 0 ) {
        return alpha*f2;
    }

    //------- cleanup data ---------------------------------------------//    
    delete[] xstencil;
    delete[] fstencil;
}
