#ifndef  SAVEPNG_H
#define  SAVEPNG_H
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
int writepng(const void * buffer, char * filename, int width, int height);
#endif
