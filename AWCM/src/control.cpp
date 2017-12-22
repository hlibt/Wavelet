#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void control(string &equation, int &max_scale, int &shift, double &threshold, int &interp_points, int &num_timesteps, double &tf, 
                double &advec_vel, double &diffusivity, string &buffer_type, int &buffer_width,
               int &buffer_height,  bool &ifwrite) {

    // read all simulation parameters from input text file
    string line;
    ifstream myfile ("input.txt");
    if ( myfile.is_open() ) {    
        while ( getline (myfile,line) ) {
            if ( line.compare(0,8,"equation") == 0 ) {
                equation = line.substr(9,string::npos);
            }
            if ( line.compare(0,9,"max_scale") == 0 ) {
                string substring = line.substr(10,string::npos);
                max_scale = stoi(substring,0);
            }
            if ( line.compare(0,5,"shift") == 0 ) {
                string substring = line.substr(6,string::npos);
                shift = stoi(substring,0);
            }
            if ( line.compare(0,9,"threshold") == 0 ) {
                string substring = line.substr(10,string::npos);
                threshold = stod(substring,0);
            }
            if ( line.compare(0,13,"interp_points") == 0 ) {
                string substring = line.substr(14,string::npos);
                interp_points = stoi(substring,0);
            }
            if ( line.compare(0,2,"tf") == 0 ) {
                string substring = line.substr(3,string::npos);
                tf = stod(substring,0);
            }
            if ( line.compare(0,13,"num_timesteps") == 0 ) {
                string substring = line.substr(14,string::npos);
                num_timesteps = stoi(substring,0);
            }
            if ( line.compare(0,9,"advec_vel") == 0 ) {
                string substring = line.substr(10,string::npos);
                advec_vel = stod(substring,0);
            }
            if ( line.compare(0,11,"diffusivity") == 0 ) {
                string substring = line.substr(12,string::npos);
                diffusivity = stod(substring,0);
            }
            if ( line.compare(0,11,"buffer_type") == 0 ) {
                buffer_type = line.substr(12,string::npos);
            }
            if ( line.compare(0,12,"buffer_width") == 0 ) {
                string substring = line.substr(13,string::npos);
                buffer_width = stoi(substring,0);
            }
            if ( line.compare(0,13,"buffer_height") == 0 ) {
                string substring = line.substr(14,string::npos);
                buffer_height = stoi(substring,0);
            }
            if ( line.compare(0,7,"ifwrite") == 0 ) {
                string substring = line.substr(8,string::npos);
                ifwrite = stoi(substring,0);
            }
        }
    } else {
        cout << "Could not open file. Check if 'input.txt' is in the correct directory." << endl; 
    }
    myfile.close();

    // create a startup message with summary of simulation parameters
    printf("==========================================================================\n");
    printf("                                                                            \n");
    printf("                        ADAPTIVE WAVELET COLLOCATION SOLVER                 \n");
    printf("                   ONE DIMENSIONAL ADVECTION-DIFFUSION EQUATION             \n");
    printf("                                                                            \n");
    printf("    Maximum wavelet level: %d \n", max_scale);
    printf("    Starting level: %d \n", shift);
    printf("    Wavelet coefficient threshold parameter: %.2e \n", threshold);
    printf("    Number of timesteps: %d \n", num_timesteps);
    printf("    Final simulation time: %3.3f \n", tf);
    printf("    Advection velocity: %3.3f \n", advec_vel);
    printf("    Diffusivity coefficient: %3.3f \n", diffusivity);
    printf("    Buffer width: %d \n", buffer_width);
    printf("    Buffer height: %d \n", buffer_height);
    printf("    Write data to files (true/false): %d \n", ifwrite);
    printf("                                                                            \n");
    printf("==========================================================================\n");

    int x;
    printf(" Enter any key and enter to continue. \n ");
    cin >> x;
    return;
}
