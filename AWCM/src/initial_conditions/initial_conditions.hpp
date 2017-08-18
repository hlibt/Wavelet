#include <math.h>
#define PI 3.14159265

class initial_condition {
    public:  
       double f(double x); 
};

// Initial condition:
double initial_condition::f(double x) {
    return cos(80.*PI*x)*exp(-64.*x*2);  
}
