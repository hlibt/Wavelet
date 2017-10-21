void RK4(double* u,double* f,bool* activPnts,double h,int n) {
    double k1, k2, k3, k4;
    for (int i=1;i<n-1;i++) {
        if (activPnts[i]==true) {
            k1=f[i];
            k2=f[i]+.5*h*k1;
            k3=f[i]+.5*h*k2;
            k4=f[i]+h*k3;
            u[i]=f[i]+(h/6.)*(k1+k2+k3+k4);
        }
    }
    u[0]=0.;
    u[n-1]=0.;
    return;
}
