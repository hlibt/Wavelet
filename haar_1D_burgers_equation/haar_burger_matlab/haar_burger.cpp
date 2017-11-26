#include <iostream>
#include <vector>
using namespace std;

main()
{
    
}

//*****************************************************************************//

void BiCGSTAB(vector<double>* Diag,vector<double>* Lowr,vector<double>* Uppr,
    vector<double>* b,vector<double>* x,double tol,int mxi)
// Purpose:  Solves the matrix equation 'Ax=b' using the Biconjugate gradient
//           stabilized method. 
//
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
{
    vector<double> r_new;
    vector<double> r_old;
    vector<double> rhat;
    vector<double> v_new;
    vector<double> v_old;
    vector<double> p_new;
    vector<double> p_old;
    vector<double> h;
    vector<double> s;
    vector<double> t;
    vector<double> x_0;
    double omega_new;
    double omega_old;
    double rho_new;
    double rho_old;
    double alpha;
    double beta;
    double error;

    // Populate initial guess array:
    for (i=1,i<=N,i++)
        x_0[0]=.0;
 
}
