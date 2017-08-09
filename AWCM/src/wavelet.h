#include <math.h>
#define PI 3.14159265

class wavelet {
    public:  

	// daughter wavelet parameters
       	double a0;
       	double aj;
      	double bjk;

	// number of recursions to determine wavelet function
	int n=6; 
	
	// define daubechies-4 filter coefficients
        double c[4] = {
         0.4829629131445341E+00, 
         0.8365163037378079E+00, 
         0.2241438680420133E+00, 
        -0.1294095225512603E+00 };

	// class function declarations
	father(double x);
	mother(double x);
	daughter(double x);	
	set_params(int j,int k);
};

double wavelet::father(double x) {
  double y;
  if ( 0 < n ) {
      y=sqrt(2.)*
        (c[0]*father(n-1,2.0*x) 
        +c[1]*father(n-1,2.0 *x-1.0) 
        +c[2]*father(n-1,2.0*x-2.0) 
        +c[3]*father(n-1,2.0*x-3.0));
  }
  else if ( 0.0 <= x && x < 1.0 ) {
    y = 1.0;
  }
  else {
    y = 0.0;
  }
  return y;
}

double wavelet::mother(double x) {
    double y;
    y=sqrt(2.0)* 
        ( pow(-1.,0.)*c[3]*father(n-1,2.0*x) 
        + pow(-1.,1.)*c[2]*father(n-1,2.0*x-1.0) 
        + pow(-1.,2.)*c[1]*father(n-1,2.0*x-2.0) 
        + pow(-1.,3.)*c[0]*father(n-1,2.0*x-3.0) );
    return y;
}

double wavelet::daughter(double x) {
    double y;
    y=pow(aj,-.5)*mother((x-bjk)/aj);    
    return y;
}

void set_params(int j,int k) {
    aj=
    
}
