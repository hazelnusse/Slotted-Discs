CC = gcc
CFLAGS = -Wall -g -O3
CXX = g++
CXXFLAGS = -Wall -g -O3

all : render simulate 

render : render.o slotted_discs.o writepng.o
	$(CXX) -lm -lgsl -lcblas -latlas -lGLU -lOSMesa -lpng -o $@ $^

simulate : simulate.o slotted_discs.o
	$(CXX) -lm -lgsl -lcblas -latlas -o $@ $^

simulate.o : simulate.cpp slotted_discs.h

slotted_discs.o : slotted_discs.cpp slotted_discs.h

render.o : render.cpp writepng.h slotted_discs.h

writepng.o : writepng.c writepng.h

clean : 
	rm -rf *.o simulate render

.PHONY: all clean
