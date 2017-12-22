void compute_field(CollocationPoint** collPnt);
void compute_derivatives( CollocationPoint** collPnt ); 
void time_integrate(CollocationPoint** collPnt,double h,string &equation, double c,double alpha);
void seed_grid(CollocationPoint** collPnt);
void wavelet_derivative(double* f,double* df,double h,int N,int derivative_order);
