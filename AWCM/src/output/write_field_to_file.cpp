#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "../CollocationPoint.hpp"
#include "../global.hpp"
using namespace std;

void write_field_to_file(CollocationPoint** collPnt,int current_timestep) {

    char u_out[45];                                                     // declare space for character array
    sprintf(u_out,"output/_soln_files/u%d.dat",current_timestep);       // append filename to include timestep number
    ofstream output;                            	
    output.open(u_out);     
    int num_active = 0;    
    CollocationPoint* tmp = new CollocationPoint[jPnts(J)];
    for (int j=0;j<=J;j++) {
        int N = jPnts(j);
        for (int i=0;i<N;i++) {  
            if ( collPnt[j][i].isMask == true ) {
                int k = indexShift(J,j,i);
                tmp[k] = collPnt[j][i];
                num_active++;
            }
        }
    }
    for (int i=0;i<jPnts(J);i++) {
        if ( tmp[i].isMask == true ) {
            output << std::fixed << std::setprecision(16) << tmp[i].x << " " << tmp[i].u << endl;    
        }
    }
    output.close(); 
    printf(" \n");
    printf("------------------------------------------ \n");
    printf("Timestep: %d \n",current_timestep);
    printf("Active points: %d out of %d \n",num_active,jPnts(J));
    printf("------------------------------------------ \n");   
    printf(" \n");
    delete[] tmp;
}
