double* BiCGSTAB(double** A,double* b,double tol,int size,int mxi);
bool chk_conv(double** A,double* y,double* b,double tolerance,int n);
double inner_product(double* a, double* b,int n);
double* ADOTX(double** A,double* x,int n);
double L2norm(double* x,int n);
