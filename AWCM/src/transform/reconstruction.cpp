#include <cmath>
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"

void reconstruction(double** x,double* solution,double** scalCoeff,double** detCoeff,int J,int numInterpPnts) {

    double** phi=new double*[J+1];
    double** psi=new double*[J+1];
    for (int j=0;j<=J;j++) {
        int N=jPnts(j);
        phi[j]=new double[N];
        psi[j]=new double[N];
    }
    for (int j=0;j<J;j++) {
        int N=jPnts(j);
        for (int l=0;l<N;l++) {
            if (j==0) scaling_subd(phi,x,j,l,J,numInterpPnts);
            detail_subd(psi,x,j,l,J,numInterpPnts);
            for (int i=0;i<jPnts(J);i++) {
                if (j==0) solution[i]+=scalCoeff[j][l]*phi[J][i];
                if (l<N-1) solution[i]+=detCoeff[j][l]*psi[J][i]; 
            }
            phi[j][l]=0.;
            psi[j][l]=0.;
        }
    }
    delete[] phi;
    delete[] psi;

    return;
}
