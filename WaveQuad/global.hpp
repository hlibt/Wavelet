extern int shift;
extern int J;
extern int interpPnts;

//------- inline function declarations ---------------------------------//              
int inline jPnts(int j) { return pow(2,j+shift) + 1; }                  // the number of grid points at level j
                                                                        //
double inline integrand( double x ) {                                   // the function to be integrated
    double f;
    f = exp( x );
/*    if ( x == 0. ) {
        f = 0.;
    } else {
        f = x * log( abs(x) );                                       // integrand
    } */
    return f;                                                           //
}                                                                       //    
