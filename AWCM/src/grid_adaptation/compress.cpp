#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void compress(CollocationPoint** collPnt, double** gridpts, double** funcpts, int** spatial_identifiers, int** level_identifiers, int& counter) {

    //--------------------------------------------------------------------------//
    // Information: This program compresses all active points in the 
    //              dyadic grid into a set of vectors. This process
    //              greatly simplifies the computation of derivatives 
    //              and time stepping on the adaptive grid. A mapping is 
    //              created to map points on the adaptive grid to the dense 
    //              vector for computations and file writing.
    //           
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              epsilon                 - threshold parameter for detail coeffs.
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    //--------------------------------------------------------------------------//     
    
    //------- initialize full vector for level J ---------------------------//
    double* xtmp = new double[jPnts(J)];                                    // vector for grid points at finest level
    double* ftmp = new double[jPnts(J)];                                    // vector for function points at finest level
    int* stmp = new int[jPnts(J)];                                          // spatial index mapping
    int* ltmp = new int[jPnts(J)];                                          // level index mapping
    bool* active = new bool[jPnts(J)];                                      // denotes active points on grid
    
    //------- represent all points at finest level -------------------------//
    counter=0;                                                              //
    for (int j=0;j<=J;j++) {                                                // loop through all other levels
        int N = jPnts(j);                                                   // number of points at level j
        for (int k=0;k<N;k++) {                                             // loop through points at level j

            //------- create vectors corresponding to finest level ---------//
            int l = indexShift(J,j,k);                                      // compute appropriate fine level index
            if ( collPnt[j][k].isMask == true ) {                           // 
                xtmp[l] = collPnt[j][k].x;                                  //
                ftmp[l] = collPnt[j][k].u;                                  //
                stmp[l] = k;                                                //
                ltmp[l] = j;                                                //
                active[l] = true;                                           //
                counter++;                                                  //
            } else {                                                        //
                active[l] = false;                                          //
            }

        }
    }

    //------- compress the grid --------------------------------------------//
    int l = 0;                                                              //
    gridpts = new double*[counter];                                          //
    funcpts = new double*[counter];                                          //
    spatial_identifiers = new int*[counter];                                 //
    level_identifiers = new int*[counter];                                   //
    for (int k=0;k<jPnts(J);k++) {                                          //
        
        //------- if points are active put them in compressed grid ---------//
        if ( active[k] == true ) {                                          //
            gridpts[l][0] = xtmp[k];                                           //
            funcpts[l][0] = ftmp[k];                                           //
            spatial_identifiers[l][0] = stmp[k];                               //
            level_identifiers[l][0] = ltmp[k];                                 //
            l++;                                                            //
        }                                                                   //
                                                                            //
    }                                                                       //

    //------- cleanup data -------------------------------------------------//
    delete[] xtmp;                                                          //
    delete[] ftmp;                                                          //
    delete[] stmp;                                                          //
    delete[] ltmp;                                                          //
    delete[] active;                                                        //

    return;
}
