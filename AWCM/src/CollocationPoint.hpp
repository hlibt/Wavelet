#include <vector>
using namespace std;

class CollocationPoint {

    //------- Public members -----------------------//------------------------------------------------------------------//
    public: 
        double x;                                   // location on the one-dimensional grid
        double u;                                   // the solution at the point
        double ux;                                  // the first derivative wrt x at the point
        double uxx;                                 // the second derivative wrt x at the point
        double scaling_coeff;                       // the scaling coefficient at the point
        double detail_coeff;                        // the wavelet coefficient at the point (if it is an odd point)
        bool isMask;                                // stores whether or not the point is in the mask at current iterate
};
