#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#define REAL float
#define MAXLEVEL 12
using namespace std;
REAL f(REAL x) { return pow(x,6.) + sin(x); };
double df(double x) { return 6.*pow(x,5.) + cos(x); }

int main() {

    // create arrays
    REAL stepsize[MAXLEVEL+1];
    REAL cntrd2pnt[MAXLEVEL+1];
    REAL cntrd4pnt[MAXLEVEL+1];
    REAL lagr_derivative[MAXLEVEL+1];
    double stencil[4];
    REAL func_pnts[4];
    double error2pnt[MAXLEVEL+1];
    double error4pnt[MAXLEVEL+1];
    double error_lagrange[MAXLEVEL+1];

    // number of grid points to test
    int numpnts = 1000;

    // compute stepsizes at all simulated wavelet levels
    for (int j=0;j<MAXLEVEL;j++) {
        stepsize[j] = pow(2.,-j-1);
    }

    // compute sample points for differentiation
    for (int i=0;i<numpnts;i++) {
        
        // compute a grid point
        REAL x = i*0.001;
        
        // compute 2-point and 4-point centered finite differences
        for (int j=0;j<MAXLEVEL;j++) {
            REAL h = stepsize[j];
            cntrd2pnt[j] = ( f(x+h) - f(x-h) ) / ( 2.*h );
            cntrd4pnt[j] = ( f(x-2.*h) - 8.*f(x-h) + 8.*f(x+h) - f(x+2.*h) ) 
                            / ( 12. * h );
        }
    
        // compute stencil points of the lagrange polynomial
        for (int j=0;j<MAXLEVEL;j++) {
            double h = pow(2.,-j-1);
            stencil[0] = i*.001 - 2.5*h;
            stencil[1] = i*0.001 + h;
            stencil[2] = i*0.001 + 2.*h;
            stencil[3] = i*0.001 + 3.*h;
            func_pnts[0] = f(stencil[0]);
            func_pnts[1] = f(stencil[1]);
            func_pnts[2] = f(stencil[2]);
            func_pnts[3] = f(stencil[3]);
            lagr_derivative[j] = 0.;
            for (int jj=0;jj<4;jj++) {
                double sum = 0.;
                for (int ii=0;ii<4;ii++) {
                    double productInner;
                    if (ii!=jj) {
                        productInner = 1.;
                        for (int m=0;m<4;m++) {
                            if (m!=ii && m!=jj) {
                                double difference = stencil[jj] - stencil[m];
                                productInner *= (i*0.001-stencil[m]) / difference;
                            }
                        }
                        double difference = stencil[jj] - stencil[ii];
                        sum += productInner / difference;
                    }
                }
                lagr_derivative[j] += func_pnts[jj] * sum;
            }
        }

        // compute errors for derivative computations
        for (int j=0;j<MAXLEVEL;j++) {
            error2pnt[j] += abs(df(x) - cntrd2pnt[j]); 
            error4pnt[j] += abs(df(x) - cntrd4pnt[j]); 
            error_lagrange[j] += abs(df(i*0.001) - lagr_derivative[j]);
        }

    }

    // compute the mean of the errors and write to file
    ofstream output;
    output.open("error.dat");
    for (int j=0;j<MAXLEVEL;j++) {
        error2pnt[j] = error2pnt[j] / numpnts;
        error4pnt[j] = error4pnt[j] / numpnts;
        error_lagrange[j] = error_lagrange[j] / numpnts;
        output << stepsize[j] << " " << error2pnt[j] << " " << error4pnt[j] << " " << error_lagrange[j] << endl;
    }
    output.close();
    return 0;
}
