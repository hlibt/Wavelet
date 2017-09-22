#include <math.h>

class initial_condition {
    public:  
       double f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
    return cos(80*M_PI*x)*exp(-64*x*x);  
//    return cos(M_PI*x)  ;
//    return x*x;
}
