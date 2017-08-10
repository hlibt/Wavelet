#include <iostream>

#include <math.h>
#include <stdio.h>
#define PI 3.14159265

class wavelet {
    public:  

	// daughter wavelet parameters
    double b0;
  	double bjk;
   	double a0;   
   	double aj;

	// class function declarations
	double father(int n,double x);
	double mother(int n,double x);
	double daughter(double x);	
	void set_params(int j,int k);

};

double wavelet::father(int n,double x) {
  
  // define daubechies-4 filter coefficients
  double c[4] = {
         0.4829629131445341E+00, 
         0.8365163037378079E+00, 
         0.2241438680420133E+00, 
        -0.1294095225512603E+00 };

  double y;
  if ( 0 < n ) {
      y=sqrt(2.)*
        (c[0]*wavelet::father(n-1,2.0*x) 
        +c[1]*wavelet::father(n-1,2.0 *x-1.0) 
        +c[2]*wavelet::father(n-1,2.0*x-2.0) 
        +c[3]*wavelet::father(n-1,2.0*x-3.0));
  }
  else if ( 0.0 <= x && x < 1.0 ) {
    y = 1.0;
  }
  else {
    y = 0.0;
  }
  return y;
}

double wavelet::mother(int n,double x) {

	// define daubechies-4 filter coefficients
    double c[4] = {
         0.4829629131445341E+00, 
         0.8365163037378079E+00, 
         0.2241438680420133E+00, 
        -0.1294095225512603E+00 };

    double y;
    y=sqrt(2.0)* 
        ( pow(-1.,0.)*c[3]*wavelet::father(n-1,2.0*x) 
        + pow(-1.,1.)*c[2]*wavelet::father(n-1,2.0*x-1.0) 
        + pow(-1.,2.)*c[1]*wavelet::father(n-1,2.0*x-2.0) 
        + pow(-1.,3.)*c[0]*wavelet::father(n-1,2.0*x-3.0) );
    return y;
}

double wavelet::daughter(double x) {
    int n=6;
    double y;
    y=pow(aj,-.5)*wavelet::mother(n,(x-bjk)/aj);    
    return y;
}

void set_params(int j,int k) {

    // largest scale
    int L=1;

    // grid boundaries
    double xl=-1.;
    double xr=1.;

	// daughter wavelet parameters
    b0=1.;
  	a0=pow(2.,-L)*(xr-xl)/b0;   
   	    
    aj=pow(2.,-(j+1.))*a0;    
    bjk=xl+aj*b0*k;   
}

