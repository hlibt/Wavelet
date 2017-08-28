#include <math.h>
#define PI 3.14159265

class initial_condition {
    public:  
       double f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
    return exp(-pow(x,2.)/pow(0.1,2));  
//    return cos(PI*x)  ;
//    return x*x;
}
