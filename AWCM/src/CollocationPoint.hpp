#include <vector>
#include "interpolation/interpolation.hpp"
using namespace std;

class CollocationPoint {

    //------- Public members -----------------------//------------------------------------------------------------------//
    public: 
        int j;                                      // the level of the collocation point
        int k;                                      // the translation parameter for the point
        double x;                                   // location on the one-dimensional grid
        double coefficient;                         // the coefficient at the point
        virtual vector<double> compute_wavelet();   // compute the scaling wavelet
        bool isMask;                                // stores whether or not the point is in the mask at current iterate
};

class ScalingPoint : public CollocationPoint {

    //------- Public members -----------------------------------------------//------------------------------------------//
    vector<double> compute_wavelet( vector<double> gridPnts , int Jmax ) {
        return scaling_subd( gridPnts, j, k, Jmax );
    };
};

class DetailPoint : public CollocationPoint {

    //------- Public members -----------------------------------------------//------------------------------------------//
    vector<double> compute_wavelet( vector<double> gridPnts , int Jmax ) {
        return detail_subd( gridPnts, j, k, Jmax );   
    };
};
