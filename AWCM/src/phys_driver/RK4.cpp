void RK4(double* u, double* f, double h, int n) {
    double k1, k2, k3, k4;
    for (int i=1;i<n;i++) {
//        k1 = f[i];
//        k2 = u[i] + .5 * h * k1;
//        k3 = u[i] + .5 * h * k2;
//        k4 = u[i] + h * k3;
//        u[i] = u[i] + (h/6.) * ( k1 + 2. * k2 + 2. * k3 + k4 );
        u[i] = u[i] + h*f[i];
    }
    u[0] = u[n-1];
    return;
}
