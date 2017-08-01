#include <math.h>
#define PI 3.14159265

class initial_condition {
    public:  
       f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
//    return sin(2.*PI*x);
    return 4.*x*(1.-x);
//    return 1.;
}
