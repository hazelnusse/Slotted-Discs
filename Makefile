all : sim render

sim : sim.o twindiscs.o
	g++ -Wall -g -lm -lgsl -lcblas -latlas -o sim sim.o twindiscs.o

sim.o : simulate.cpp twindiscs.h
	g++ -Wall -g -c simulate.cpp -o sim.o

twindiscs.o : twindiscs.cpp twindiscs.h
	g++ -Wall -g -c twindiscs.cpp

render : render.o twindiscs.o writepng.o
	g++ -Wall -g -lm -lgsl -lcblas -latlas -lGLU -lOSMesa -lpng -o render render.o twindiscs.o writepng.o

render.o : render.cpp writepng.h slotted_discs.h
	g++ -Wall -g -c render.cpp

writepng.o : writepng.cpp writepng.h
	g++ -Wall -g -c writepng.cpp

clean : 
	rm -rf *.o simulate render sim
