double lagrInterp_old(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN);
double lagrInterp(double xeval, CollocationPoint** collPnt, int j, int i, int numPnts);
double lagrInterpD1(double x, CollocationPoint** collPnt, int j, int i, int numPnts);
double lagrInterpD2(double x,double* gridPnts,double* funcPnts,int i,int n,int maxN);
void scaling_subd(double** f,double** x,int j,int m,int Jmax,int N);
void detail_subd(double** f,double** x,int j,int m,int Jmax,int N);
double kronecker_delta(int k, int m);
