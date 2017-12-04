#include <vector>
using namespace std;

class CollocationPoint {

    //------- Public members -----------------------//------------------------------------------------------------------//
    public: 
        int j;                                      // the level of the collocation point
        int k;                                      // the translation parameter for the point
        double x;                                   // location on the one-dimensional grid
        double scaling_coeff;                       // the scaling coefficient at the point
        double detail_coeff;                        // the wavelet coefficient at the point (if it is an odd point)
        bool isMask;                                // stores whether or not the point is in the mask at current iterate
};
