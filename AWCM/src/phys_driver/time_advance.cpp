#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "../CollocationPoint.hpp"
#include "../grid_adaptation/grid_adaptation.hpp"
#include "../global.hpp"
#include "phys_driver.hpp"

void time_advance(CollocationPoint** collPnt,double h,string equation, double c,double alpha,
                        string boundary_type, double left_bc, double right_bc, int nactive) {

    //------- initialize compression vectors -------------------------------//
    double* gridPts = new double[nactive];                                  // vector of active grid points
    double* funcPts = new double[nactive];                                  // vector of functional values at those points
    int* spatialMapPts = new int[nactive];                                  // integer array of spatial mapping values
    int* levelMapPts = new int[nactive];                                    // integer array of level mapping values

    //------- compress the grid --------------------------------------------//
    compress(collPnt,gridPts,funcPts,spatialMapPts,levelMapPts,nactive);    // compress grid and output number of active points

    //------- Advance interior points in time ------------------------------//
    RK2(gridPts,funcPts,nactive,h,alpha,c,equation);                        // move solution variables forward

    //------- Compute boundary values --------------------------------------//
    if ( boundary_type.compare(0,9,"derichlet") == 0 ) {                    //
        collPnt[0][0].u = left_bc;                                          //
        collPnt[0][jPnts(0)-1].u = right_bc;                                //
    }                                                                       //

    //------- decompress the grid ------------------------------------------//
    decompress(collPnt,funcPts,spatialMapPts,levelMapPts,nactive);          // go from compressed grid back to dyadic-adaptive

    //------- cleanup data -------------------------------------------------//
    delete[] gridPts;                                                       //
    delete[] funcPts;                                                       //
    delete[] spatialMapPts;                                                 //
    delete[] levelMapPts;                                                   //

    return;
}
