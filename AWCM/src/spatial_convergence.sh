#!/bin/bash

#--------- Information --------------------------------------------#
# This Bash script automates the generation of input 
# files (from a template) for your PDE solver, runs the code,
# and stores output in desired locations in an organized fashion. 
# Performing convergence studies can be a pain in the neck (literally) 
# as the code must be rerun many dozens of times with different 
# parameters. The less efficient method to do this is to modify the 
# source code itself, essentially putting the entire code in a loop
# while grid spacing or time steps are varied. But this is not
# the most convenient solution as the code is not modular, and 
# the parameters to be changed can be in various places. This code
# puts everything in one convenient place and runs as many options
# as you like. Please modify and use as you please. Many upgrades
# to the code are possible and more sophisticated methods have
# been devised, but are not as easy to find and use. 
# The names of my input files are hardcoded as '_in_32pts.txt' for 
# 32 cells for example but this could easily be modified. The names 
# of the folders for the solution data are named like '_run_32pts/'.
# Note that on Mac, the 'sed' command needs a backup file specified.
# The C++ PDE solver I use takes in input file and output folder 
# names as 'argv' from the command line.
#
# Developed by Brandon Gusto. 
#------------------------------------------------------------------#



#---------- Functions ---------------------------------------------#

Create () {
if [ ! -d $1 ]; then
    mkdir $1;
fi
}

#------------------------------------------------------------------#



#---------- PDE code options and parameters -----------------------#

# Reconstruction options recognized by the control file:
#... "upwind_difference" or
#... "central_difference" or
#... "muscl2" or
#... "linear_upwind_difference" or
#... "piecewise_parabolic" or
#... "mono_piecewise_parabolic"

# Declare methods that you want to run:
declare -a METHDS=( "upwind_difference"
                    "central_difference"
                    "muscl2"
                    "piecewise_parabolic"
                    "mono_piecewise_parabolic" )

# Timestep options recognized by the control file:
#... "explicit" if you want to explicitly declare a timestep, or
#... "cfl" if you want timestep controlled by CFL number

# Timestep option:
TSO="explicit"

# Courant number (ignored if "explicit" declared):
CFL=0.5

# Timestep size (ignored if "cfl" declared):
DT=0.001

# Number of quadrature nodes per cell:
NQ=2

# Final simulation time:
TF=10.0
#6.28318530717958647692528677

# Constant advection velocity:
VEL=1.0

# Left and right ends of the domain (2*pi is 6.28318530717958647692528677):
XL=0.0
XR=1.0

# Write-frequency:
WF=1

#------------------------------------------------------------------#



#---------- Convergence study parameters --------------------------#

# Specify minimum number of grid points to start:
NMIN=8

# Factor to increase number of grid cells by each run:
FACTOR=2

# Number of runs:
NRUNS=6

# Folder for input files:
INDIR=input
Create ${INDIR}

# Folder for output files:
OUTDIR=output
Create ${OUTDIR}

#------------------------------------------------------------------#



#---------- Perform the convergence study -------------------------#

# Loop through the methods:
for m in ${METHDS[@]}; do

    # Create the input file and run the program:
    let NPTS=${NMIN}
    for ((i=1;i<=NRUNS;i++)) do

        # Display the current method and number of cells:
        echo "Currently running '${m}' with $NPTS points..."

        # Create folder for input files for each method
        DIR1=${INDIR}/${m}
        Create ${DIR1}

        # Create input file name with desired parameters:
        FILE=${INDIR}/${m}/_in_${NPTS}pts.txt

	    # Copy the template input file:
        cp ${INDIR}/input_template.txt ${FILE}

        # Edit the file parameters:
        sed -i .bak "s/^timestep_select.*/timestep_select ${TSO}/" ${FILE}
        sed -i .bak "s/^cfl_number.*/cfl_number ${CFL}/" ${FILE}
        sed -i .bak "s/^step_size.*/step_size ${DT}/" ${FILE}
        sed -i .bak "s/^number_cells.*/number_cells ${NPTS}/" ${FILE}
        sed -i .bak "s/^number_quadrature_points.*/number_quadrature_points ${NQ}/" ${FILE}
        sed -i .bak "s/^final_time.*/final_time ${TF}/" ${FILE}
        sed -i .bak "s/^advection_speed.*/advection_speed ${VEL}/" ${FILE}
        sed -i .bak "s/^left_grid_boundary.*/left_grid_boundary ${XL}/" ${FILE}
        sed -i .bak "s/^right_grid_boundary.*/right_grid_boundary ${XR}/" ${FILE}
        sed -i .bak "s/^reconstruction_scheme.*/reconstruction_scheme ${m}/" ${FILE}
        sed -i .bak "s/^write_frequency.*/write_frequency ${WF}/" ${FILE}
        
        # Create folders for individual run output data:
        DIR2=${OUTDIR}/${m}
        Create ${DIR2}

        DIR3=${DIR2}/_run_${NPTS}pts
        Create ${DIR3}

        # Run the executable (ppa) and pass input file and output folder names:
        ./ppa ${FILE} ${DIR3}

        # Increment the number of grid cells:
        let NPTS=${NPTS}*${FACTOR}

    done

done
