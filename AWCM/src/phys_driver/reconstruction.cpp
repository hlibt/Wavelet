#include <cmath>
#include <vector>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
using namespace std;

void reconstruction(CollocationPoint** collPnt) {

    //--------------------------------------------------------------------------//
    // Information: This program uses the active wavelets to reconstruct the 
    //              original function.  
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              J                       - maximum grid level (global)
    //--------------------------------------------------------------------------//     

    double** gridPnts = new double*[J+1];                                       // create arrays for grid points to pass through scaling_subd
    double** phi = new double*[J+1];                                            // arrays for the scaling function
    double** psi = new double*[J+1];                                            // arrays for the detail wavelet
    for (int j=0;j<=J;j++) {                                                    // create columns
        int N = jPnts(j);                                                       //
        gridPnts[j] = new double[N];                                            //
        phi[j] = new double[N];                                                 //
        psi[j] = new double[N];                                                 //
    }                                                                           //
    for (int j=0;j<=J;j++) {                                                    //
        int N = ( jPnts(j) - 1 ) / 2;                                           //
        for (int k=-N;k<=N;k++) {                                      	        //
            gridPnts[j][k+N] = 2. *  pow(2.,-(j+shift)) * k;                    // x-locations of each collocation point
        }                                                       	            //
    }                                                                           //

    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int l=0;l<N;l++) {
            collPnt[j][l].u = 0.;
        }
    }
                                                                                //
    for (int j=0;j<J;j++) {                                                     // 
        int N = jPnts(j);                                                         
        for (int l=0;l<N;l++) {
            if (j==0) scaling_subd(phi,gridPnts,j,l,J,interpPnts);              // generate the scaling function if at level j=0
            if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {              // if at odd point and the point is point is in the mask,
                detail_subd(psi,gridPnts,j,l,J,interpPnts);                     // then generate the detail wavelet at that point
            }                                                                   //
            for (int i=0;i<jPnts(J);i++) {
                if (j==0) {
                    collPnt[J][i].u += collPnt[j][l].scaling_coeff * phi[J][i];
                }
                if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {
                    collPnt[J][i].u += collPnt[j+1][2*l+1].detail_coeff * psi[J][i]; 
                }
            }
        }
    }
    delete[] phi;
    delete[] psi;
    delete[] gridPnts;

    return;
}

