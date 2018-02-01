void seed_grid(CollocationPoint** collPnt);
void compute_field(CollocationPoint** collPnt);
void time_advance(CollocationPoint** collPnt,double h,string equation, double c,
     double alpha,string boundary_type, double left_bc, double right_bc);
void RK2(double* gridpts, double* funcpts, int nactive, double h, double alp, double c, string equation);
double rhs(double* gridpts, double* funcpts, int nactive, int i, double mod, double alpha, double c, string equation);
