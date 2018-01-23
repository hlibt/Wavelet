#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#define REAL double
#define MAXLEVEL 12
using namespace std;
typedef complex<REAL> dcomp;
REAL f(REAL x) { return pow(x,3.) - cos(x); };
REAL df(REAL x) { return 3.*pow(x,2.) + sin(x); }

int main() {

    // create arrays
    REAL stepsize[MAXLEVEL+1];
    REAL cntrd2pnt[MAXLEVEL+1];
    REAL cntrd4pnt[MAXLEVEL+1];
    double re2pnt[MAXLEVEL+1];
    double error2pnt[MAXLEVEL+1];
    double error4pnt[MAXLEVEL+1];
    double error_re[MAXLEVEL+1];

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

        // compute richardson extrapolation of the 2-point finite difference
        for (int j=0;j<MAXLEVEL;j++) {
            REAL h = stepsize[j];
            REAL a = ( f(x+h) - f(x-h) ) / (2.*h);
            REAL b = ( f(x+2.*h) - f(x-2.*h) ) / (4.*h);
            double y = pow(a,0.5);
            double z = pow(b,0.5);
            re2pnt[j] = ( 3.*a + (y-z)*(y+z) ) / 3.; 
        }

        // compute errors for derivative computations
        for (int j=0;j<MAXLEVEL;j++) {
            error2pnt[j] += abs(df(x) - cntrd2pnt[j]); 
            error4pnt[j] += abs(df(x) - cntrd4pnt[j]); 
            error_re[j] += abs(df(x) - re2pnt[j]);
        }

    }

    // compute the mean of the errors and write to file
    ofstream output;
    output.open("RE.dat");
    for (int j=0;j<MAXLEVEL;j++) {
        error2pnt[j] = error2pnt[j] / numpnts;
        error4pnt[j] = error4pnt[j] / numpnts;
        error_re[j] = error_re[j] / numpnts;
        output << stepsize[j] << " " << error2pnt[j] << " " << error4pnt[j] << " " << error_re[j] << endl;
    }
    output.close();
    return 0;
}
