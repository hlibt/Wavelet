#include <math.h>

class initial_condition {
    public:  
       double f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
    return cos(80*M_PI*x)*exp(-64*x*x);  
//    double x0=1./3.;
//    double v=pow(10.,-2.);
//    return -tanh((x+x0)/(2.*v))+exp(-pow(64.,2.)*pow(x-x0,2.));
//    return sin(M_PI*x)  ;
}
