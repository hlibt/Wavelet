#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
#include "../interpolation/interpolation.hpp"

double rhs(CollocationPoint** collPnt, int j, int i, double mod, string equation) { 

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

    //------- Modify the gridpoint's stencil values --------------------//
    double* gridPnts = new double[2*interpPnts];                        // the interpolating stencil
    double* funcPnts = new double[2*interpPnts];                        // the function evaluated at stencil points
    double f0 = collPnt[j][i].u + mod;                                  // solution at time subinterval
    int p = (i-1)/2 - interpPnts + 1;                                   // starting grid index
    for (int k=0;k<2*interpPnts;k++) {                                  // loop through stencil
        gridPnts[k] = collPnt[j-1][p].x;                                //
        funcPnts[k] = collPnt[j-1][p].u + mod;                          //
        p++;                                                            //
    }

    //------- Compute first-order spatial derivative -------------------//
    double f1 = 0.;                                                     // the first derivative of u
    double x = collPnt[j][i].x;                                         // evaluation point
    for (p=0;p<2*interpPnts;p++) {
        double sum = 0.;
        for (int k=0;k<2*interpPnts;k++) {
            if (k!=p) {
                double productUpper = 1.;
                for (int l=0;l<2*interpPnts;l++) {
                    if ( l!=k && l!=p ) { 
                        productUpper *= ( x - gridPnts[l] );
                    }
                }
                sum += productUpper;
            }
        }
        double productLower = 1.;
        for (int k=0;k<2*interpPnts;k++) {
            if (k!=p) { 
                productLower *= ( gridPnts[p] - gridPnts[k] );
            }
        }
        f1 += funcPnts[p] * sum / productLower;
    }

    //------- Compute second-order spatial derivative ------------------//
    double f2 = 0.;                                                     // the second derivative of u
    for (p=0;p<2*interpPnts;p++) {                                      //
        double sumOuter = 0.;                                           //
        for (int l=0;l<2*interpPnts;l++) {                              //
            if (l!=i) {                                                 //
                double sumInner = 0.;                                   //
                for (int m=0;m<2*interpPnts;m++) {                      //
                    if (m!=i && m!=l) {                                 //
                        double product=1.;                              //
                        for (int k=0;k<2*interpPnts;k++) {              //
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
