all: awcm

awcm: awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o lagrInterpD1.o
	g++ -o awcm awcm.o scaling_subd.o detail_subd.o fwd_trans.o lagrInterp.o lagrInterpD1.o -I.

awcm.o: awcm.cpp initial_conditions/initial_conditions.hpp
	g++ -c awcm.cpp -I.

scaling_subd.o: interpolation/scaling_subd.cpp 
	g++ -c interpolation/scaling_subd.cpp -I. -I/interpolation/

detail_subd.o: interpolation/detail_subd.cpp 
	g++ -c interpolation/detail_subd.cpp -I -I/interpolation/

fwd_trans.o: transform/fwd_trans.cpp
	g++ -c transform/fwd_trans.cpp -I -I/transform

lagrInterp.o: interpolation/lagrInterp.cpp
	g++ -c interpolation/lagrInterp.cpp -I -I/interpolation

lagrInterpD1.o: interpolation/lagrInterpD1.cpp
	g++ -c interpolation/lagrInterpD1.cpp -I -I/interpolation   

clean:
	rm *o

