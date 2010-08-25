simulate : simulate.o slotted_discs.o
	gcc -Wall -g -o simulate simulate.o slotted_discs.o

simulate.o : simulate.c 
	gcc -Wall -g -c simulate.c

slotted_discs.o : slotted_discs.c slotted_discs.h
	gcc -Wall -g -c slotted_discs.c

