extern int shift;
extern int J;
extern int interpPnts;

//------- inline function declarations ---------------------------------//              
int inline jPnts(int j) { return pow(2,j+shift) + 1; }                  // the number of grid points at level j
                                                                        //
double inline integrand( double x ) {                                   // the function to be integrated
//    double f = log( 1 - cos( x ) );                                     // singluar integrand
    double f = x*x;
    return f;                                                           //
}                                                                       //    
