all : simulate render sim

simulate : simulate.o slotted_discs.o
	gcc -Wall -g -lm -lgsl -lcblas -latlas -o simulate simulate.o slotted_discs.o

sim : sim.o twindiscs.o
	g++ -Wall -g -lm -lgsl -lcblas -latlas -o sim sim.o twindiscs.o

sim.o : simulate.cpp twindiscs.h
	g++ -Wall -g -c simulate.cpp -o sim.o

twindiscs.o : twindiscs.cpp twindiscs.h
	g++ -Wall -g -c twindiscs.cpp

simulate.o : simulate.c
	gcc -Wall -g -c simulate.c

slotted_discs.o : slotted_discs.c slotted_discs.h
	gcc -Wall -g -c slotted_discs.c

render : render.o slotted_discs.o writepng.o
	gcc -Wall -g -lm -lgsl -lcblas -latlas -lGLU -lOSMesa -lpng -o render render.o slotted_discs.o writepng.o

render.o : render.c writepng.h slotted_discs.h
	gcc -Wall -g -c render.c

writepng.o : writepng.c writepng.h
	gcc -Wall -g -c writepng.c

clean : 
	rm -rf *.o simulate render sim
