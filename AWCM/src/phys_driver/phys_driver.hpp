int thresholding(CollocationPoint** collPnt, double epsilon);
void reconstruction(CollocationPoint** collPnt);
void compute_derivatives( CollocationPoint** collPnt ); 
void compute_rhs(double* rhs, double* ux, double c, int active);
void RK4(double* u, double* f, double h, int n);
