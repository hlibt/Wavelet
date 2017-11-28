#include <cmath>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"

void reconstruction(double* solution,CollocationPoint collPnt) {

    double** phi = new double*[J+1];                                                // arrays for scaling function phi
    double** psi = new double*[J+1];                                                // arrays for wavelets psi
    for (int j=0;j<=J;j++) {                                                        // create the columns
        int N = jPnts(j);                                                           //
        phi[j] = new double[N];                                                     //
        psi[j] = new double[N];                                                     //
    }                                                                               //
    for (int j=0;j<J;j++) {                                                         // loop through all scales up to the finest
        int N = jPnts(j);                                                           // number of points at level j
        for (int l=0;l<N;l++) {                                                     // loop through all points at level j
            if (j==0) scaling_subd(phi,x,j,l,J,interpPnts);                         // compute the scaling functions
            if (mask[j+1][2*l+1]==true && l<N-1) detail_subd(psi,x,j,l,J,numInterpPnts);
//            if (l<N-1) detail_subd(psi,x,j,l,J,numInterpPnts);
            for (int i=0;i<jPnts(J);i++) {
                if (j==0) {
                    solution[i]+=scalCoeff[j][l]*phi[J][i];
                }
                if (l<N-1) {
                    solution[i]+=detCoeff[j][l]*psi[J][i]; 
                }
                phi[j][i]=0.;
                psi[j][i]=0.;
            }
        }
    }
    delete[] phi;
    delete[] psi;

    return;
}
