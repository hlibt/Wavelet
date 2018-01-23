#include <cmath>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
#include "phys_driver.hpp"
using namespace std;

void compute_field(CollocationPoint** collPnt) {

    //--------------------------------------------------------------------------//
    // Information: compute_field.cpp reconstructs the field by
    //              summing the active wavelets and wavelet coefficients.
    //              Derivatives of the function are computed by differentiating
    //              the active wavelets. 
    //
    // Input: 
    //              collPnt                 - the matrix of collocation point objects 
    //                                        at all levels. Composed of other variables
    //              collPnt.scaling_coeff   - the scaling coefficients
    //              collPnt.detail_coeff    - the detail coefficients
    //              epsilon                 - threshold parameter for detail coeffs.
    //              J                       - maximum grid level (global variable)
    //              interpPnts              - half the number of interp. points (global)
    // Output:
    //              collPnt.u               - the field u(x)
    //              collPnt.ux              - the first derivative of u(x)
    //              collPnt.uxx             - the second derivative of u(x)
    //--------------------------------------------------------------------------//     

    //------- Allocate space for necessary arrays ------------------------------//
    int N = jPnts(J);                                                           // number of points at maximum level of resolution 
    double* phi = new double[N];                                                // arrays for the scaling function
    double* psi = new double[N];                                                // arrays for the detail wavelet
    double* u_reconstructed = new double[N];                                    // reconstructed solution variable

    //------- Compute interval length for spatial derivatives ------------------//
    double h = abs( collPnt[J][1].x - collPnt[J][0].x );                        //

    //------- Set to zero all reconstructed variables --------------------------//
    for (int i=0;i<jPnts(J);i++) {                                              //
        u_reconstructed[i] = 0.;                                                //
    }                                                                           //

    //------- Loop through all levels collecting significant wavelets ----------//
    for (int j=0;j<J;j++) {                                                     // 
        int N = jPnts(j);                                                       //
        for (int l=0;l<N;l++) {                                                 // 
            if (j==0) {                                                         // j=0 for the points corresponding to scaling functions
                scaling_subd(collPnt,phi,j,l);                                  // generate the scaling function if at level j=0
            }
            if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {              // if at odd point and the point is point is in the mask,
                detail_subd(collPnt,psi,j,l);                                   // then generate the detail wavelet at that point
            }                                                                   //
            for (int i=0;i<jPnts(J);i++) {
                if (j==0) {
                    u_reconstructed[i] += collPnt[j][l].scaling_coeff * phi[i];
                }
                if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {
                    u_reconstructed[i] += collPnt[j+1][2*l+1].detail_coeff * psi[i]; 
                }
            }
        }
    }

    //------- Transfer results to the adaptive dyadic grid ---------------------//
    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int l=0;l<N;l++) {
            if ( collPnt[j][l].isMask == true ) {
                int k = indexShift(J,j,l);
                collPnt[j][l].u = u_reconstructed[k];
            }
        }
    }
    
    delete[] phi;
    delete[] psi;
    delete[] u_reconstructed;

    return;
}

