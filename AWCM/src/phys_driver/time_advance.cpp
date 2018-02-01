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

<<<<<<< HEAD
    //------- Advance interior points in time --------------------------------------//
    for (int j=0;j<=J;j++) {                                                        // search all levels for points in the mask
        int N = jPnts(j);                                                           // number of points at level j
        for (int i=0;i<N;i++) {                                                     // loop through points at level j
 
            //------- Advance point if it is in mask -------------------------------//
            if ( collPnt[j][i].isMask == true ) {                                   // determine if point is in mask
                RK2(collPnt,j,i,h,equation);                                        // call time integration scheme
            }                                                                       //
=======
    //------- initialize compression vectors -------------------------------//
    double* gridPts;                                                        // vector of active grid points
    double* funcPts;                                                        // vector of functional values at those points
    int* spatialMapPts;                                                     // integer array of spatial mapping values
    int* levelMapPts;                                                       // integer array of level mapping values
    int na;                                                                 //

    //------- compress the grid --------------------------------------------//
    compress(collPnt,gridPts,funcPts,spatialMapPts,levelMapPts,na);         // compress grid and output number of active points
>>>>>>> 2f6e487d56af862323ff001ca16a8666c02e1e39

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
