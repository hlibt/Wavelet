#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "../CollocationPoint.hpp"
#include "../grid_adaptation/grid_adaptation.hpp"
#include "../global.hpp"
#include "phys_driver.hpp"

void time_advance(CollocationPoint** collPnt,double h,string equation, double c,double alpha,
                        string boundary_type, double left_bc, double right_bc) {

    //------- initialize compression vectors -------------------------------//
    double* gridPts;                                                        // vector of active grid points
    double* funcPts;                                                        // vector of functional values at those points
    int* spatialMapPts;                                                     // integer array of spatial mapping values
    int* levelMapPts;                                                       // integer array of level mapping values
    int na;                                                                 //

    //------- compress the grid --------------------------------------------//
    compress(collPnt,gridPts,funcPts,spatialMapPts,levelMapPts,na);         // compress grid and output number of active points

    //------- Advance interior points in time ------------------------------//
    RK2(gridPts,funcPts,na,h,alpha,c,equation);                             // move solution variables forward

    //------- decompress the grid ------------------------------------------//
    decompress(collPnt,funcPts,spatialMapPts,levelMapPts,na);               // go from compressed grid back to dyadic

    //------- Compute boundary values --------------------------------------//
    for (int j=0;j<=J;j++) {                                                //
        if ( boundary_type.compare(0,9,"derichlet") == 0 ) {                //
            collPnt[j][0].u = left_bc;                                      //
            collPnt[j][jPnts(j)-1].u = right_bc;                            //
        }                                                                   //
    }                                                                       //

    //------- cleanup data -------------------------------------------------//
    delete[] gridPts;                                                       //
    delete[] funcPts;                                                       //
    delete[] spatialMapPts;                                                 //
    delete[] levelMapPts;                                                   //

    return;
}
