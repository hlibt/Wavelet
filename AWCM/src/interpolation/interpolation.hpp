double lagrInterp(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN);
double lagrInterpD1(double x,double* gridPnts,double* funcPnts,int n);
void scaling_subd(double** f,double** x,int j,int m,int Jmax,int N);
void detail_subd(double** f,double** x,int j,int m,int Jmax,int N);
double lagrInterp(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN);
double lagrInterpD1(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN);
double kronecker_delta(int k, int m);
