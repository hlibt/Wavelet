void RK4(double* f,double h,int n) {
    double k1, k2, k3, k4;
    for (int i=0;i<n;i++) {
        k1=f[i];
        k2=f[i]+.5*h*k1;
        k3=f[i]+.5*h*k2;
        k4=f[i]+h*k3;
        f[i]=f[i]+(h/6.)*(k1+k2+k3+k4);
    }
    return;
}
