void compute_field(CollocationPoint** collPnt);
double rhs(CollocationPoint** collPnt, int j, int i, double mod, string equation);
void time_integrate(CollocationPoint** collPnt,double h,string equation, double c,double alpha,
                        string boundary_type, double left_bc, double right_bc);
void seed_grid(CollocationPoint** collPnt);
void RK2(CollocationPoint** collPnt, int j, int i, double h, string equation);
void wavelet_derivative(double* f,double* df,double h,int N,int derivative_order);
