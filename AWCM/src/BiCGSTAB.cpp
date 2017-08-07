double* BiCGSTAB(double** A,double* b,double tol,int size,int mxi) {
    // Purpose:  Solves the matrix equation 'Ax=b' using
    //          the Biconjugate gradient stabilized method. 
    // Modified: July 03, 2017
    //
    // Author: Brandon Gusto
    //
    // Parameters: Input --> Diag, the elements on the diagonal of the matrix 'A'
    //             Input --> Lowr, the elements on the lower diagonal of the matrix 'A'
    //             Input --> Uppr, the elements on the upper diagonal of the matrix 'A'
    //             Input --> b, the elements in the vector 'b' in the RHS of 'Ax=b'
    //             Input --> x, the solution vector
    //             Input --> tol, the tolerance with which the BiCGSTAB method solves the system
    //             Input --> mxi, the maximum number of iterations
    //             Local -->
    double* r_new=new double[size];
    double* r_old=new double[size];
    double* rhat=new double[size];
    double* v_new=new double[size];    
    double* v_old=new double[size];
    double* p_new=new double[size];
    double* p_old=new double[size];
    double* h=new double[size];
    double* s=new double[size];
    double* t=new double[size];
    double* x0=new double[size];
    double* tmp=new double[size];
    double omega_new;
    double omega_old;
    double rho_new;
    double rho_old;
    double alpha;
    double beta;

    for (int i=0;i<size;i++) x0[i]=1.0;             // populate initial solution guess
    tmp=ADOTX(A,x0,size);                           // calculate A*x0;
    for (int i=0;i<size;i++) r_old[i]=b[i]-tmp[i];  // start the residual with initial guess
    for (int i=0;i<size;i++) rhat[i]=r_old[i];      // choose arbitrary vector rhat
    rho_old=1.; omega_old=1.; alpha=1.;             // initial parameters
    for (int i=0;i<size;i++) {
        v_old[i]=0.;
        p_old[i]=0.;
    }
    int k=1;
    bool iterate=true;
    do {
        rho_new=inner_product(rhat,r_old,size);
        beta=(rho_new/rho_old)*(alpha/omega_old);
        for (int i=0;i<size;i++) p_new[i]=r_old[i]+beta*(p_old[i]-omega_old*v_old[i]);
        v_new=ADOTX(A,p_new,size);
        alpha=rho_new/inner_product(rhat,v_new,size);
        for (int i=0;i<size;i++) h[i]=x0[i]+alpha*p_new[i];
        if (chk_conv(A,h,b,tol,size)==true) {
            iterate=false;
            for (int i=0;i<size;i++) x0[i]=h[i];
            cout << "Check 1 \n";
        }
        else {
            for (int i=0;i<size;i++) s[i]=r_old[i]-alpha*v_new[i];
            t=ADOTX(A,s,size);
            omega_new=inner_product(t,s,size)/inner_product(t,t,size);
            for (int i=0;i<size;i++) x0[i]=h[i]+omega_new*s[i];
            if (chk_conv(A,x0,b,tol,size)==true) {
                iterate=false;
                cout << "Check 2 \n";
            }
            else {
                for (int i=0;i<size;i++) r_new[i]=s[i]-omega_new*t[i];
            }
        }
        for (int i=0;i<size;i++) {
            r_old[i]=r_new[i];
            p_old[i]=p_new[i];
            v_old[i]=v_new[i];
        }            
        rho_old=rho_new;
        omega_old=omega_new;
        k++;
        if (k>mxi) {
            iterate=false;
            cout << "Max iterations reached \n";
        }
    }while(iterate==true);
    return x0;
}
