CC = gcc
CFLAGS = -Wall -g -O3
CXX = g++
CXXFLAGS = -Wall -g -O3

all : render simulate 

render : render.o twindiscs.o writepng.o
	$(CXX) -lm -lgsl -lcblas -latlas -lGLU -lOSMesa -lpng -o $@ $^

simulate : simulate.o twindiscs.o
	$(CXX) -lm -lgsl -lcblas -latlas -o $@ $^

simulate.o : simulate.cpp twindiscs.h

twindiscs.o : twindiscs.cpp twindiscs.h

render.o : render.cpp writepng.h twindiscs.h

writepng.o : writepng.c writepng.h

clean : 
	rm -rf *.o simulate render

.PHONY: all clean
