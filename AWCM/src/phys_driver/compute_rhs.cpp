
void compute_rhs(double* rhs, double* ux, double c, int active) {

    for (int i=0;i<active;i++) {
            rhs[i] = -c * ux[i];
    }
    
}
