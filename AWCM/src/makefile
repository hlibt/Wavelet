all: run

run: awcm.o wavelet_generation/scaling_subd.o
	g++ awcm.o wavelet_generation/scaling_subd.o -o run

awcm.o: awcm.cpp
	g++ awcm.cpp

wavelet_generation/scalind_subd.o: wavelet_generation/scaling_subd.cpp
	g++ wavelet_generation/scaling_subd.cpp

clean:
	rm *o run

