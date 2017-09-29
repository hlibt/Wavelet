#include <math.h>

class initial_condition {
    public:  
       double f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
//    return cos(80*M_PI*x)*exp(-64*x*x);  
    return -tanh((x+1./3.)/(2.*pow(10.,-2.)))+exp(-pow(64.,2.)*pow(x-1./3.,2.));
//    return cos(M_PI*x)  ;
//    return x*x;
}
