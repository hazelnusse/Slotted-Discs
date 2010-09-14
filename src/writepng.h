#ifdef __cplusplus
extern "C" {
#endif

#ifndef WRITEPNG_H 
#define WRITEPNG_H

#include <png.h>
#include <stdio.h>
#include <stdlib.h>
int writepng(const void * buffer, char * filename, int width, int height);

#endif

#ifdef __cplusplus
}
#endif
