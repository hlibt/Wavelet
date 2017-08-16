#include <math.h>
#define PI 3.14159265

class initial_condition {
    public:  
       double f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
    return 4.*x*(1.-x);
}
