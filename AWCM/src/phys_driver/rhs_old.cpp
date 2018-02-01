#include <cmath>
#include <iostream>
#include <string>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
#include "../interpolation/interpolation.hpp"
using namespace std;

double rhs(CollocationPoint** collPnt, int j, int i, double mod, double alpha, double c, string equation) { 

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

    //------- initialize vectors ---------------------------------------//
    int N = 2*interpPnts;                                               // number of points in the stencil
    int leftPnt;                                                        // leftmost index of stencil
    int rightPnt;                                                       // rightmost index of stencil
    int p;                                                              // leftmost grid index of the stencil
    double* gridPnts = new double[N];                                   // the interpolating stencil
    double* funcPnts = new double[N];                                   // the function evaluated at stencil points
    double x = collPnt[j][i].x;                                         // evaluation point of derivative
    double f0 = collPnt[j][i].u + mod;                                  // solution at time subinterval

    //------- find adequate stencil for differentiation ----------------//
    while(1) {                                                          // break loop when adequate stencil found
        if ( j!=J ) {                                                   // only check point above if not at finest level
            if ( i>0 && collPnt[j+1][2*i+1].isMask == false ) {         // check if point above is not in mask
                leftPnt = i - interpPnts + 1;                           // starting index of stencil to use
                rightPnt = interpPnts + i;                              // end index of stencil
                break;                                                  //  
            } else if ( i<jPnts(j)-1 &&                                 //
                    collPnt[j+1][2*i-1].isMask == false ) {             // check adjacent point
                leftPnt = i - interpPnts;                               // starting index of stencil to use
                rightPnt = interpPnts + i - 1;                          // end index of stencil
                break;                                                  //
            } else {                                                    // go up a level
                j++;                                                    // search one grid level up
                i*=2;                                                   // adjust spatial index accordingly
            }                                                           //
        } else {                                                        // 
            leftPnt = i - interpPnts + 1;                               //
            rightPnt = interpPnts + i;                                  //
        }                                                               //
    }                                                                   //
    while ( leftPnt < 0 ) {                                             //
        leftPnt++;                                                      //
        rightPnt++;                                                     //
    }                                                                   //
    while ( rightPnt > jPnts(j)-1 ) {                                   //
        leftPnt--;                                                      //
        rightPnt--;                                                     //
    }                                                                   //

    //------- modify the gridpoint's stencil values --------------------//
    for (int k=0;k<N;k++) {                                             // loop through stencil
        if ( collPnt[j][leftPnt+k].isMask == true ) {                   //
            gridPnts[k] = collPnt[j][leftPnt+k].x;                      //
            funcPnts[k] = collPnt[j][leftPnt+k].u + mod;                // add modification to solution variable
        } else {
            cout << "Stencil Error." << endl;
        }
    }

    //------- compute first-order spatial derivative -------------------//
    double f1 = 0.;                                                     // the first derivative of u
    for (p=0;p<N;p++) {
        double sum = 0.;
        for (int k=0;k<N;k++) {
            if (k!=p) {
                double productUpper = 1.;
                for (int l=0;l<N;l++) {
                    if ( l!=k && l!=p ) { 
                        productUpper *= ( x - gridPnts[l] );
                    }
                }
                sum += productUpper;
            }
        }
        double productLower = 1.;
        for (int k=0;k<N;k++) {
            if (k!=p) { 
                productLower *= ( gridPnts[p] - gridPnts[k] );
            }
        }
        f1 += funcPnts[p] * sum / productLower;
    }

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
                                product *= ( x - gridPnts[k] ) /        //
                                          ( gridPnts[i] - gridPnts[k]); //
                            }                                           //
                        }                                               //
                        sumInner += product / (gridPnts[i]-gridPnts[m]);//
                    }                                                   //
                }                                                       //
                sumOuter += sumInner / (gridPnts[i]-gridPnts[l]);       //
            }                                                           //
        }                                                               //
        f2 += funcPnts[i] * sumOuter;                                   //
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
    
}
