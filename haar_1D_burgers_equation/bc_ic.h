#include <math.h>
#define PI 3.14159265

class BCIC {
    public:
    
    // Function declarations:
    double f(double x);         // the initial condition
    double fx(double x);        // first derivative of initial condition
    double fxx(double x);       // second derivative of initial condition
    double f1(double t);        // boundary condition on the left end
    double f1t(double t);       // derivative of boundary function d( f1(t) ) / dt
    double f2(double t);        // boundary condition on the right end
    double f2t(double t);       // derivative of boundary function d( f2(t) ) / dt
};

// Initial condition:
double BCIC::f(double x) {
    return sin(2.*PI*x);
}

// First derivative of initial function:
double BCIC::fx(double x) {
    return 2.*PI*cos(2.*PI*x);
}

// Second derivative of initial function:
double BCIC::fxx(double x) {
    return -4.*pow(PI,2.)*sin(2.*PI*x);
}

// Left boundary condition:
double BCIC::f1(double t) {
    return 0.;
}

// Left boundary condition derivative:
double BCIC::f1t(double t) {
    return 0.;
}

// Right boundary condition:
    double BCIC::f2(double t) {
    return 0.;
}

// Right boundary condition derivative:
    double BCIC::f2t(double t) {
    return 0.;
}

