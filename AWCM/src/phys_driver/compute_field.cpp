#include <cmath>
#include <vector>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"
using namespace std;
void diffWave(double** f,double** df,double h,int J,int N,int derivative_order);

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
    double** gridPnts = new double*[J+1];                                       // create arrays for grid points to pass through scaling_subd
    double** phi = new double*[J+1];                                            // arrays for the scaling function
    double** psi = new double*[J+1];                                            // arrays for the detail wavelet
    double** Dphi=new double*[J+1];
    double** Dpsi=new double*[J+1];
    double** DDphi=new double*[J+1];
    double** DDpsi=new double*[J+1];
    double* u_reconstructed = new double[jPnts(J)];                             // reconstructed solution variable
    double* ux_reconstructed = new double[jPnts(J)];                            //
    double* uxx_reconstructed = new double[jPnts(J)];                           //
    for (int j=0;j<=J;j++) {                                                    // create columns
        int N = jPnts(j);                                                       //
        gridPnts[j] = new double[N];                                            //
        phi[j] = new double[N];                                                 //
        psi[j] = new double[N];                                                 //
        Dphi[j] = new double[N];                                                //
        Dpsi[j] = new double[N];                                                //
        DDphi[j] = new double[N];                                               //
        DDpsi[j] = new double[N];                                               //
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
        uxx_reconstructed[i] = 0.;
    }

    for (int j=0;j<J;j++) {                                                     // 
        int N = jPnts(j);                                                         
        for (int l=0;l<N;l++) {
            if (j==0) {
                scaling_subd(phi,gridPnts,j,l,J,interpPnts);              // generate the scaling function if at level j=0
                diffWave(phi,Dphi,abs(gridPnts[J][1]-gridPnts[J][0]),J,jPnts(J),1);   //
                diffWave(phi,DDphi,abs(gridPnts[J][1]-gridPnts[J][0]),J,jPnts(J),2);   //
            }
            if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {              // if at odd point and the point is point is in the mask,
                detail_subd(psi,gridPnts,j,l,J,interpPnts);                     // then generate the detail wavelet at that point
                diffWave(psi,Dpsi,abs(gridPnts[J][1]-gridPnts[J][0]),J,jPnts(J),1);   //
                diffWave(psi,DDpsi,abs(gridPnts[J][1]-gridPnts[J][0]),J,jPnts(J),2);   //
            }                                                                   //
            for (int i=0;i<jPnts(J);i++) {
                if (j==0) {
                    u_reconstructed[i] += collPnt[j][l].scaling_coeff * phi[J][i];
                    ux_reconstructed[i] += collPnt[j][l].scaling_coeff * Dphi[J][i];
                    uxx_reconstructed[i] += collPnt[j][l].scaling_coeff * DDphi[J][i];
                }
                if ( collPnt[j+1][2*l+1].isMask == true && l < N-1 ) {
                    u_reconstructed[i] += collPnt[j+1][2*l+1].detail_coeff * psi[J][i]; 
                    ux_reconstructed[i] += collPnt[j+1][2*l+1].detail_coeff * Dpsi[J][i];
                    uxx_reconstructed[i] += collPnt[j+1][2*l+1].detail_coeff * DDpsi[J][i];
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
                collPnt[j][l].uxx = uxx_reconstructed[k];
            }
        }
    }
    
    delete[] phi;
    delete[] psi;
    delete[] Dphi;
    delete[] Dpsi;
    delete[] DDphi;
    delete[] DDpsi;
    delete[] u_reconstructed;
    delete[] ux_reconstructed;
    delete[] uxx_reconstructed;
    delete[] gridPnts;

    return;
}

void diffWave(double** f,double** df,double h,int J,int N,int derivative_order) {
    if ( derivative_order == 1 ) {
        df[J][0]=(-3.*f[J][0]+4.*f[J][1]-f[J][2])/(2.*h);
        for (int i=1;i<N-1;i++) {
            df[J][i]=(f[J][i+1]-f[J][i-1])/(2.*h);
        }
        df[J][N-1]=(3.*f[J][N-1]-4.*f[J][N-2]+f[J][N-2])/(2.*h);
    }
    if ( derivative_order == 2 ) {
        df[J][0] = ( -f[J][3] + 4.*f[J][2] - 5.*f[J][1] + 2.*f[J][0] ) / ( h*h );
        for (int i=1;i<N-1;i++) {
            df[J][i] = ( f[J][i+1] - 2.*f[J][i] + f[J][i-1] ) / ( h*h );
        }
        df[J][N-1] = ( -f[J][N-4] + 4.*f[J][N-3] - 5.*f[J][N-2] + 2.*f[J][N-1] ) / ( h*h );
    }
    return;
}
