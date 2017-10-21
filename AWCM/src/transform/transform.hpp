void fwd_trans(double** x,double** scalCoeff,double** detCoeff,int Jmax,int N);
void reconstruction(double** x,double* solution,double** scalCoeff,double** detCoeff,bool** mask,int J,int numInterpPnts);
