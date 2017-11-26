all: quadrature

quadrature: quadrature.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o lagrInterp_old.o thresholding.o wavelet_integrals.o integrate.o
	g++ -o quadrature quadrature.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o lagrInterp_old.o thresholding.o wavelet_integrals.o integrate.o -I.

quadrature.o: quadrature.cpp global.hpp CollocationPoint.hpp phys_driver/phys_driver.hpp transform/transform.hpp interpolation/interpolation.hpp
	g++ -c quadrature.cpp -I.

scaling_subd.o: interpolation/scaling_subd.cpp
	g++ -c interpolation/scaling_subd.cpp -I. -I/phys_driver

detail_subd.o: interpolation/detail_subd.cpp
	g++ -c interpolation/detail_subd.cpp -I. -I/phys_driver

fwd_trans.o: transform/fwd_trans.cpp
	g++ -c transform/fwd_trans.cpp -I. -I/transform/

lagrInterp.o: interpolation/lagrInterp.cpp
	g++ -c interpolation/lagrInterp.cpp -I. -I/interpolation/

lagrInterp_old.o: interpolation/lagrInterp_old.cpp
	g++ -c interpolation/lagrInterp_old.cpp -I. -I/interpolation/

thresholding.o: phys_driver/thresholding.cpp
	g++ -c phys_driver/thresholding.cpp -I. -I/phys_driver/

wavelet_integrals.o: phys_driver/wavelet_integrals.cpp
	g++ -c phys_driver/wavelet_integrals.cpp -I. -I/phys_driver/

integrate.o: phys_driver/integrate.cpp
	g++ -c phys_driver/integrate.cpp -I. -I/phys_driver/

clean:
	rm *o


