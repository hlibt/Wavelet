all: awcm

awcm: awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o
	g++ -o awcm awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o -I.

awcm.o: awcm.cpp initial_conditions/initial_conditions.hpp
	g++ -c awcm.cpp -I.

scaling_subd.o: wavelet_generation/scaling_subd.cpp wavelet_generation/wavelet_generation.hpp
	g++ -c wavelet_generation/scaling_subd.cpp

detail_subd.o: wavelet_generation/detail_subd.cpp wavelet_generation/wavelet_generation.hpp
	g++ -c wavelet_generation/detail_subd.cpp

fwd_trans.o: transform/fwd_trans.cpp
	g++ -c transform/fwd_trans.cpp

lagrInterp.o: interpolation/lagrInterp.cpp
	g++ -c interpolation/lagrInterp.cpp

clean:
	rm *o awcm

