#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double lagrInterp(double x, double* gridPts, double* funcPts, int i, int npts, int Nlagr);

double fct(double x) {
    //return -tanh((x+1./3.)/(2.*pow(10.,-2.)))+exp(-pow(64.,2.)*pow(x-1./3.,2.));
    return cos(M_PI*x);
}

int main()
{
	int npts = 65;                                  // npts := number of points of true function
	int neval = 64;
	double* f = new double [npts];                  // f := true function to be approximated
	double* x = new double [npts];                  // x := grid of size npts for f to be evaluated at
	double* xeval = new double [neval];
	double xmin = 0.;
    double xmax = 1.;
	double h = (xmax - xmin) / (npts-1);
	double heval = (xmax - xmin) / (neval-1);
	int Nlagr = 4;                                  // Nlagr := number of points in the interpolation stencil

	double pi = acos(-1.);

	for (int i=0; i < npts; i++) {
		x[i] = i*h;
		f[i] = fct(x[i]);
	}

    // evaluate polynomials at construction points
    for (int i=0; i < npts; i++) {
		double ff = lagrInterp(x[i], x, f, i, npts, Nlagr);
//		printf("gridPt[%d]= %14.7f , funcPts[%d]= %14.7f, lagr[%d]= %14.7f\n", i, x[i], i, f[i], i, ff);
	}

	// Evaluate function at neval pts and compute the error norm
	double norm = 0.;
	for (int i=0; i < npts-1; i++) {
        xeval[i] = .5*(x[i]+x[i+1]);  // redefine xeval to be the midpoint between x[i]'s (just like in the actual awcm)
		double feval = lagrInterp(xeval[i], x, f, i, npts, Nlagr);
		double err = feval - fct(xeval[i]); // determine error 
		norm += err * err;
	}
	norm = sqrt(norm); // L2 norm
	printf("norm= %21.14f\n", norm);
    return 0;
}


// Given grid points gridPts[0] ... gridPts[n-1]
// and a point x, create a Lagrange interpolation scheme with N points
double lagrInterp(double x, double* gridPts, double* funcPts, int i, int npts, int Nlagr)
{
// Establish Lagrange interpolant centered at gridPts[i]. If centering is not possible shift
// as needed to make it possible. If Nlagr is odd, gridPts[i] should be the center. 
// If Nlagr is even, which points should the interpolant use (Brian can implement)

/*    if (Nlagr % 2 == 0) {
		// Not implemented yet
		return -1.;
	} */

	// find the point gridPts[i] closest to x
	// NOT DONE

	int leftPt  = i - (Nlagr-1)/2;
	int rightPt = i + (Nlagr-1)/2 + 1; // When looping ,go from leftPt to rightPt  (c-style)
//	printf("leftPt,rightPt= %d, %d\n", leftPt, rightPt);

	if (leftPt < 0) {
		rightPt -= leftPt;
		leftPt   = 0;
	} else if (rightPt > (npts-1)) {
		leftPt += (npts-1-rightPt);
		rightPt = npts-1;
	}
	printf("leftPt,rightPt= %d, %d\n", leftPt, rightPt);
   
    double sum = 0.;

	// Cost is O(N^2). Isn't that very expensive?
    for (int l = leftPt; l <= rightPt; l++) {
        double product = 1.;
        for (int k = leftPt; k <= rightPt; k++) {
            if (k != l) {
                product *= (x-gridPts[k]) / (gridPts[l]-gridPts[k]);
            }
        }
        sum+=product*funcPts[l];
    }
    return sum;
}
