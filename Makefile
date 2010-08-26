all : simulate

simulate : simulate.o slotted_discs.o
	gcc -Wall -g -lm -lgsl -lcblas -latlas -o simulate simulate.o slotted_discs.o

simulate.o : simulate.c 
	gcc -Wall -g -c simulate.c

slotted_discs.o : slotted_discs.c slotted_discs.h
	gcc -Wall -g -c slotted_discs.c

clean : 
	rm -rf *.o simulate
