#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void write2file(CollocationPoint** collPnt,int timestep) {

    //------- Output the field -----------------------------------------//
    char u_out[45];                                                     // declare space for character array
    sprintf(u_out,"output/_soln_files/u%d.dat",timestep);               // append filename to include timestep number
    ofstream output;                            	
    output.open(u_out);     
    int num_active = 0;    
    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {  
            if ( collPnt[j][i].isMask == true ) {
                output << std::fixed << std::setprecision(16) <<        //
                    collPnt[j][i].x << " " << collPnt[j][i].u << endl;  //    
                if ( collPnt[j][i].isOdd == true ) num_active++;        // only count active wavelets
            }
        }
    }
    output.close(); 

    //------- Output the magnitude of the detail coefficients to file --//
    char coeffs_out[45];                                                // declare space for character array
    sprintf(coeffs_out,"output/_coeff_files/coeff%d.dat",timestep);     // append filename to include timestep number
    output.open(coeffs_out);     
    for (int j=1;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {  
            if ( collPnt[j][i].isOdd == true && collPnt[j][i].isMask == true ) {
                output << std::fixed << std::setprecision(16) << collPnt[j][i].x 
                        << " " << j << " " << abs( collPnt[j][i].detail_coeff ) << endl;    
                num_active++;
            }
        }
    }
    output.close(); 

    //------- Print information to the screen --------------------------//
    printf(" \n");
    printf("------------------------------------------ \n");
    printf("Timestep: %d \n",timestep);
    printf("Active points: %d out of %d \n",num_active,jPnts(J));
    printf("------------------------------------------ \n");   
    printf(" \n");
}
