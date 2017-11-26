extern int shift;
extern int J;
extern int interpPnts;

// Inline functions:
int inline jPnts(int j) { return pow(2,j+shift) + 1; }
int inline indexShift(int jstar, int j, int k); 
