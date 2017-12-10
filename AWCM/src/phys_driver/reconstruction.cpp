#include <cmath>
#include <vector>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
using namespace std;
void diffWave(double** f,double** df,double h,int J,int N);

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
    double** Dphi=new double*[J+1];
    double** Dpsi=new double*[J+1];
    double* u_reconstructed = new double[jPnts(J)];                             // reconstructed solution variable
    double* ux_reconstructed = new double[jPnts(J)];                            //
    for (int j=0;j<=J;j++) {                                                    // create columns
        int N = jPnts(j);                                                       //
        gridPnts[j] = new double[N];                                            //
        phi[j] = new double[N];                                                 //
        psi[j] = new double[N];                                                 //
        Dphi[j] = new double[N];                                                //
        Dpsi[j] = new double[N];                                                //
    }                                                                           //
    for (int j=0;j<=J;j++) {                                                    //
        int N = ( jPnts(j) - 1 ) / 2;                                           //
        for (int k=-N;k<=N;k++) {                                      	        //
            gridPnts[j][k+N] = 2. *  pow(2.,-(j+shift)) * k;                    // x-locations of each collocation point
        }                                                       	            //
    }                                                                           //

    for (int i=0;i<jPnts(J);i++) {
        u_reconstructed[i] = 0.;
        ux_reconstructed[i] = 0.;
    }

    for (int j=0;j<J;j++) {                                                     // 
        int N = jPnts(j);                                                         
        for (int l=0;l<N;l++) {
            if (j==0) {
                scaling_subd(phi,gridPnts,j,l,J,interpPnts);              // generate the scaling function if at level j=0
                diffWave(phi,Dphi,abs(gridPnts[J][1]-gridPnts[J][0]),J,jPnts(J));   //
            }
            if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {              // if at odd point and the point is point is in the mask,
                detail_subd(psi,gridPnts,j,l,J,interpPnts);                     // then generate the detail wavelet at that point
                diffWave(psi,Dpsi,abs(gridPnts[J][1]-gridPnts[J][0]),J,jPnts(J));   //
            }                                                                   //
            for (int i=0;i<jPnts(J);i++) {
                if (j==0) {
                    u_reconstructed[i] += collPnt[j][l].scaling_coeff * phi[J][i];
                    ux_reconstructed[i] += collPnt[j][l].scaling_coeff * Dphi[J][i];
                }
                if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {
                    u_reconstructed[i] += collPnt[j+1][2*l+1].detail_coeff * psi[J][i]; 
                    ux_reconstructed[i] += collPnt[j+1][2*l+1].detail_coeff * Dpsi[J][i];
                }
            }
        }
    }

    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int l=0;l<N;l++) {
            if ( collPnt[j][l].isMask == true ) {
                int k = indexShift(J,j,l);
                collPnt[j][l].u = u_reconstructed[k];
                collPnt[j][l].ux = ux_reconstructed[k];
            }
        }
    }
    
    delete[] phi;
    delete[] psi;
    delete[] Dphi;
    delete[] Dpsi;
    delete[] u_reconstructed;
    delete[] gridPnts;

    return;
}

void diffWave(double** f,double** df,double h,int J,int N) {
    df[J][0]=(-3.*f[J][0]+4.*f[J][1]-f[J][2])/(2.*h);
    for (int i=1;i<N-1;i++) {
        df[J][i]=(f[J][i+1]-f[J][i-1])/(2.*h);
    }
    df[J][N-1]=(3.*f[J][N-1]-4.*f[J][N-2]+f[J][N-2])/(2.*h);
    return;
}
