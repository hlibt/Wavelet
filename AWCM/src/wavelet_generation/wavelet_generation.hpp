double* scaling_subd(double** x,int j,int m,int Jmax,int N);
double* detail_subd(double** x,int j,int m,int Jmax,int N);
double lagrange_interp(double eval_point,double* x,int i,int N1,int N2);
double kronecker_delta(int k, int m);
double neville(double y,double* x,double* f,int n);

