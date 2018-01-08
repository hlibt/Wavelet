CXX = g++
DEBUG_FLAG =

all: awcm

objects = awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o \
			compute_field.o thresholding.o time_integrate.o	\
			seed_grid.o adjacent_zone.o control.o reconstruction_check.o \
			write2file.o wavelet_derivative.o

awcm: $(objects) 
	${CXX} -o awcm $(objects)

awcm.o: awcm.cpp conditions/conditions.hpp phys_driver/phys_driver.hpp transform/transform.hpp CollocationPoint.hpp global.hpp
	${CXX} ${DEBUG_FLAG} -c awcm.cpp -I.

scaling_subd.o: interpolation/scaling_subd.cpp 
	${CXX} ${DEBUG_FLAG} -c interpolation/scaling_subd.cpp -I. -I/interpolation/

detail_subd.o: interpolation/detail_subd.cpp 
	${CXX} ${DEBUG_FLAG} -c interpolation/detail_subd.cpp -I. -I/interpolation/

lagrInterp.o: interpolation/lagrInterp.cpp
	${CXX} ${DEBUG_FLAG} -c interpolation/lagrInterp.cpp -I. -I/interpolation/

lagrInterp_old.o: interpolation/lagrInterp_old.cpp
	${CXX} ${DEBUG_FLAG} -c interpolation/lagrInterp_old.cpp -I. -I/interpolation/

lagrInterpD1.o: interpolation/lagrInterpD1.cpp
	${CXX} ${DEBUG_FLAG} -c interpolation/lagrInterpD1.cpp -I. -I/interpolation/

lagrInterpD2.o: interpolation/lagrInterpD2.cpp
	${CXX} ${DEBUG_FLAG} -c interpolation/lagrInterpD2.cpp -I. -I/interpolation/

fwd_trans.o: transform/fwd_trans.cpp
	${CXX} ${DEBUG_FLAG} -c transform/fwd_trans.cpp -I. -I/transform/

thresholding.o: grid_adaptation/thresholding.cpp
	${CXX} ${DEBUG_FLAG} -c grid_adaptation/thresholding.cpp -I. -I/phys_driver/

adjacent_zone.o: grid_adaptation/adjacent_zone.cpp
	${CXX} ${DEBUG_FLAG} -c grid_adaptation/adjacent_zone.cpp

reconstruction_check.o: grid_adaptation/reconstruction_check.cpp
	${CXX} ${DEBUG_FLAG} -c grid_adaptation/reconstruction_check.cpp

seed_grid.o: phys_driver/seed_grid.cpp
	${CXX} ${DEBUG_FLAG} -c phys_driver/seed_grid.cpp -I. -I/phys_driver

compute_field.o: phys_driver/compute_field.cpp
	${CXX} ${DEBUG_FLAG} -c phys_driver/compute_field.cpp -I. -I/phys_driver/

wavelet_derivative.o: phys_driver/wavelet_derivative.cpp
	${CXX} ${DEBUG_FLAG} -c phys_driver/wavelet_derivative.cpp -I. -I/phys_driver/

time_integrate.o: phys_driver/time_integrate.cpp
	${CXX} ${DEBUG_FLAG} -c phys_driver/time_integrate.cpp -I. -I/phys_driver/

control.o: control.cpp
	${CXX} ${DEBUG_FLAG} -c control.cpp -I.

write2file.o: output/write2file.cpp
	${CXX} ${DEBUG_FLAG} -c output/write2file.cpp -I. -I/output/

clean:
	rm awcm $(objects)
