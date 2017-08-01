#include <math.h>

class wavelet {
    public: 
        
        // Member class function declarations:
        double phi(double x); 
        double psi
        double phi_x(double x);             // first derivative of scaling function phi  
        double psi_x(double x);             // first derivative of detail function psi
        void set_params(int a, int b);
};

// Set the parameters for each new detail wavelet:
void Haar::set_params(int a, int b) {
    kk=static_cast<double>(a);
    mm=static_cast<double>(b);
    Xi1=kk/mm;
    Xi2=(kk+.5)/mm;
    Xi3=(kk+1.)/mm;   
}
