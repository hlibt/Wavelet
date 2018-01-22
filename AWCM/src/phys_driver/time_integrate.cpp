#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void time_advance(CollocationPoint** collPnt,double h,string equation, double c,double alpha,
                        string boundary_type, double left_bc, double right_bc) {

    //------- Advance interior points in time --------------------------------------//
    for (int j=1;j<=J;j++) {
        int N = jPnts(j);
        for (int i=1;i<N-1;i++) {

            //------- Advance point if it is in mask -------------------------------//
            if ( collPnt[j][i].isMask == true ) {
                RK2(collPnt,j,i,h,equation);
            }

        }






{
                if ( equation.compare(0,7,"burgers") == 0 ) {
                    collPnt[j][i].u = ( - collPnt[j][i].u * collPnt[j][i].ux + 
                               alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
                } else if ( equation.compare(0,16,"modified_burgers") == 0 ) {
                    collPnt[j][i].u = ( - ( collPnt[j][i].u + c ) * collPnt[j][i].ux 
                            + alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
                } else if ( equation.compare(0,19,"advection_diffusion") == 0 ) {
                    collPnt[j][i].u = ( - c * collPnt[j][i].ux +  
                                alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
                } else if ( equation.compare(0,9,"advection") == 0 ) {
                    collPnt[j][i].u = ( - c * collPnt[j][i].ux ) * h + 
                                        collPnt[j][i].u;
                } else if ( equation.compare(0,9,"diffusion") == 0 ) {
                    collPnt[j][i].u = alpha * collPnt[j][i].uxx * h + collPnt[j][i].u;
                }
            }
        }   
        if ( boundary_type.compare(0,9,"derichlet") == 0 ) {
            collPnt[j][0].u = left_bc;
            collPnt[j][N-1].u = right_bc;
        } else if ( boundary_type.compare(0,8,"periodic") == 0 ) {
            int i = 0;
            if ( equation.compare(0,7,"burgers") == 0 ) {
                collPnt[j][N-1].u = ( - collPnt[j][N-1].u * collPnt[j][N-1].ux + 
                           alpha * collPnt[j][N-1].uxx ) * h + collPnt[j][N-1].u;
            } else if ( equation.compare(0,16,"modified_burgers") == 0 ) {
                collPnt[j][i].u = ( - ( collPnt[j][N-1].u + c ) * collPnt[j][i].ux 
                        + alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
            } else if ( equation.compare(0,19,"advection_diffusion") == 0 ) {
                collPnt[j][i].u = ( - c * collPnt[j][i].ux +  
                            alpha * collPnt[j][i].uxx ) * h + collPnt[j][i].u;
            } else if ( equation.compare(0,9,"advection") == 0 ) {
                collPnt[j][N-1].u = ( - c * collPnt[j][N-1].ux ) * h + 
                                    collPnt[j][N-1].u;
            } else if ( equation.compare(0,9,"diffusion") == 0 ) {
                collPnt[j][i].u = alpha * collPnt[j][i].uxx * h + collPnt[j][i].u;
            }
            collPnt[j][0].u = collPnt[j][N-1].u;
        }
    }
    return;
}
