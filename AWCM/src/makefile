all: awcm

awcm: awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o lagrInterp_old.o lagrInterpD1.o lagrInterpD2.o reconstruction.o thresholding.o RK4.o
	g++ -o awcm awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o lagrInterp_old.o lagrInterpD1.o lagrInterpD2.o reconstruction.o thresholding.o RK4.o -I.

awcm.o: awcm.cpp conditions/conditions.hpp phys_driver/phys_driver.hpp transform/transform.hpp CollocationPoint.hpp
	g++ -c awcm.cpp -I.

scaling_subd.o: interpolation/scaling_subd.cpp 
	g++ -c interpolation/scaling_subd.cpp -I. -I/interpolation/

detail_subd.o: interpolation/detail_subd.cpp 
	g++ -c interpolation/detail_subd.cpp -I. -I/interpolation/

fwd_trans.o: transform/fwd_trans.cpp
	g++ -c transform/fwd_trans.cpp -I. -I/transform/

lagrInterp.o: interpolation/lagrInterp.cpp
	g++ -c interpolation/lagrInterp.cpp -I. -I/interpolation/

lagrInterp_old.o: interpolation/lagrInterp_old.cpp
	g++ -c interpolation/lagrInterp_old.cpp -I. -I/interpolation/

lagrInterpD1.o: interpolation/lagrInterpD1.cpp
	g++ -c interpolation/lagrInterpD1.cpp -I. -I/interpolation/

lagrInterpD2.o: interpolation/lagrInterpD2.cpp
	g++ -c interpolation/lagrInterpD2.cpp -I. -I/interpolation/

thresholding.o: phys_driver/thresholding.cpp
	g++ -c phys_driver/thresholding.cpp -I. -I/phys_driver/

reconstruction.o: phys_driver/reconstruction.cpp
	g++ -c phys_driver/reconstruction.cpp -I. -I/phys_driver/

RK4.o: phys_driver/RK4.cpp
	g++ -c phys_driver/RK4.cpp -I. -I/phys_driver/

clean:
	rm *o


