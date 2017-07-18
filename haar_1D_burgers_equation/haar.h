#include <math.h>

class Haar {
    public:

        // Parameters:
        int M;              // number of collocation points
        int J;              // maximum level of resolution
        int i;              // wavelet coefficient index
        int j;              // indicates level of wavelet
        double k;           // translation parameter
        double m;           // 
        double Xi1; 
        double Xi2;
        double Xi3;

        // Member class function declarations:
        double h(double x); 
        double p(double x);     // p(x) =  integral h(x) dx
        double q(double x);     // q(x) = integral p(x) dx
        void set_params(int a, int b);
        double h1(double x);
        double p1(double x);
        double q1(double x);
        double q_tilda(int i);
};

// Member function definitions:
double Haar::h(double x) {
    if (x>=Xi1 && x<Xi2)
        return 1.;
    else if (x>=Xi2 && x<Xi3)
        return -1.;
    else 
        return 0.;
}

// first integral of the haar function:
double Haar::p(double x) {
    if (x>=Xi1 && x<Xi2)
        return x-Xi1;
    else if (x>=Xi2 && x<Xi3)
        return Xi3-x;
    else 
        return 0.;
}

// second integral of the haar function:
double Haar::q(double x) {
    if (x>=0 && x<Xi1)
        return 0.;
    else if (x>=Xi1 && x<Xi2)
        return pow(.5*(x-Xi1),2.);
    else if (x>=Xi2 && x<Xi3)
        return 1./(4.*pow(m,2.))-pow(.5*(Xi3-x),2.);
    else if (x>=Xi3 && x<=1)
        return 1./(4.*pow(m,2.));
}

// set the parameters for each new detail wavelet:
void Haar::set_params(int a, int b) {
    k=static_cast<double>(a);
    m=static_cast<double>(b);
    Xi1=k/m;
    Xi2=(k+.5)/m;
    Xi3=(k+1.)/m;   
    i=a+b+1;
}

// scaling function:
double Haar::h1(double x) {
    if (x>=0 && x<1)
        return 1.;
    else
        return 0.;
}
 
// first integral of scaling function:
double Haar::p1(double x) {
    if (x>=0 && x<1)
        return x;
    else
        return 0.;
}

// second integral of scaling function:
double Haar::q1(double x) {
    if (x>=0 && x<1)
        return .5*pow(x,2.);
    else
        return 0.;
}

// q function evaluated at boundary x=1:
double Haar::q_tilda(int i) {
    if (i==1)
        return .5;
    else if (i>1)
        return 1./(4.*pow(m,2.));
}
