extern int shift;
extern int J;
extern int interpPnts;

//------- inline function declarations ---------------------------------//              
int inline jPnts(int j) { return pow(2,j+shift) + 1; }                  // the number of grid points at level j
                                                                        //
int inline indexShift(int jstar, int j, int k) {                        // represents index k at level j, at the desired level jstar
        int gridMultplr = pow(2,jstar-j);                               // multiplier parameter to get k to level jstar
        return k * gridMultplr;                                         // output what k would be at finest level jstar
}                                                                       //
                                                                        //
double inline init_condition(double x) {                                //
    double f = 0.;                                                      //
//    f = - sin( M_PI * x ); 
    f = -tanh( (x+.5) / (2*.01) );
    return f;                                                           //
}                                                                       //

