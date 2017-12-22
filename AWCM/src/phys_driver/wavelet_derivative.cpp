#include <cmath>
#include <vector>
#include "../CollocationPoint.hpp"
#include "../interpolation/interpolation.hpp"
#include "../global.hpp"

void wavelet_derivative(double* f,double* df,double h,int N,int derivative_order) {

    //------- Compute first derivative of wavelet --------------------------------------//
    if ( derivative_order == 1 ) {
        for (int i=0;i<N-1;i++) {
            if ( i == 0 ) {
                df[i] = ( -147*f[i] + 360*f[i+1] - 450*f[i+2] + 400*f[i+3] - 225*f[i+4]
                        + 72*f[i+5] - 10*f[i+6] ) / (60*h);
            } else if ( i == 1 ) {
                df[i] = ( -10*f[i-1] - 77*f[i] + 150*f[i+1] - 100*f[i+2] + 50*f[i+3]
                            - 15*f[i+4] + 2*f[i+5] ) / (60*h);
            } else if ( i == 2 ) {
                df[i] = ( 2*f[i-2] - 24*f[i-1] - 35*f[i] + 80*f[i+1] - 
                                30*f[i+2] + 8*f[i+3] - f[i+4] ) / (60.*h);
            } else if ( i == N - 1 ) {
                df[i] = ( 10*f[i-6] - 72*f[i-5] + 225*f[i-4] - 400*f[i-3] + 450*f[i-2] - 
                            360*f[i-1] + 147*f[i] ) / (60*h);
            } else if ( i == N - 2 ) {
                df[i] = ( -2*f[i-5] + 15*f[i-4] - 50*f[i-3] + 100*f[i-2] - 150*f[i-1]
                        + 77*f[i] + 10*f[i+1] ) / (60*h); 
            } else if ( i == N - 3 ) {
                df[i] = ( f[i-4] - 8*f[i-3] + 30*f[i-2] - 80*f[i-1]
                        + 35*f[i] + 24*f[i+1] - 2*f[i+2] ) / (60*h); 
            } else {
                df[i] = ( 45.*( f[i+1] - f[i-1] ) - 9.*( f[i+2] - f[i-2] ) + 
                        ( f[i+3] - f[i-3] ) ) / (60.*h);
            }
        }
    }

    //------- Compute second derivative of wavelet ------------------------------------//
    if ( derivative_order == 2 ) {
        for (int i=0;i<N-1;i++) {
            if ( i == 0 ) {
                df[i] = (812*f[i+0]-3132*f[i+1]+5265*f[i+2]-5080*f[i+3]+2970*f[i+4]-
                        972*f[i+5]+137*f[i+6])/(180*1.0*h*h);
            } else if ( i == 1 ) {
                df[i] = (137*f[i-1]-147*f[i+0]-255*f[i+1]+470*f[i+2]-285*f[i+3]+
                        93*f[i+4]-13*f[i+5])/(180*1.0*h*h);
            } else if ( i == 2 ) {
                df[i] = (-13*f[i-2]+228*f[i-1]-420*f[i+0]+200*f[i+1]+15*f[i+2]-
                        12*f[i+3]+2*f[i+4])/(180*1.0*h*h);
            } else if ( i == N - 1 ) {
                df[i] = (2*f[i-4]-12*f[i-3]+15*f[i-2]+200*f[i-1]-420*f[i+0]+
                        228*f[i+1]-13*f[i+2])/(180*1.0*h*h);
            } else if ( i == N - 2 ) {
                df[i] = (-13*f[i-5]+93*f[i-4]-285*f[i-3]+470*f[i-2]-255*f[i-1]-
                        147*f[i+0]+137*f[i+1])/(180*1.0*h*h);
            } else if ( i == N - 3 ) {
                df[i] = (137*f[i-6]-972*f[i-5]+2970*f[i-4]-5080*f[i-3]+5265*f[i-2]-
                        3132*f[i-1]+812*f[i+0])/(180*1.0*h*h); 
            } else {
                df[i] = (2*f[i-3]-27*f[i-2]+270*f[i-1]-490*f[i+0]+270*f[i+1]-27*f[i+2]+
                        2*f[i+3])/(180*1.0*h*h); 
            }
        }
    }
    return;
}
