#include <cmath>
#include <vector>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"

void wavelet_derivative(double* f,double* df,double h,int N,int derivative_order) {

    if ( derivative_order == 1 ) {
        df[0] = (-3.*f[0]+4.*f[1]-f[2]) / (2.*h);
        for (int i=1;i<N-1;i++) {
            df[i] = (f[i+1]-f[i-1]) / (2.*h);
        }
        df[N-1] = (3.*f[N-1]-4.*f[N-2]+f[N-2]) / (2.*h);
    }
    if ( derivative_order == 2 ) {
        df[0] = ( -f[3] + 4.*f[2] - 5.*f[1] + 2.*f[0] ) / ( h*h );
        for (int i=1;i<N-1;i++) {
            df[i] = ( f[i+1] - 2.*f[i] + f[i-1] ) / ( h*h );
        }
        df[N-1] = ( -f[N-4] + 4.*f[N-3] - 5.*f[N-2] + 2.*f[N-1] ) / ( h*h );
    }

    return;
}
