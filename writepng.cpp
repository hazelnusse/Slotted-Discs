/*
 * =====================================================================================
 *
 *       Filename:  writepng.c
 *
 *    Description:  Grabs pixels from OpenGL buffer and saves them to RGBA32 png
 *
 *         Author:  Dale Lukas Peterson
 *        Company:  University of California Davis
 *
 * =====================================================================================
 */

#include "writepng.h"

int writepng(const void * buffer, char * filename, int width, int height)
{
  // Open a file to write the png for
  FILE * fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "capture: Couln't open output file \"%s\"", filename);
    return 1;
  }

  // Initialize PNG write structure
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                NULL, NULL, NULL); 
  if (!png_ptr) {
    fprintf(stderr, "capture: Can't initialize png_ptr");
    return 1;
  }
  // Initialize PNG info pointer
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
     png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
     fprintf(stderr, "capture: Can't initialze info_ptr");
     return 1;
  }
  // Initialize PNG error jump
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    fprintf(stderr, "capture: Unknown error");
    return 1;
  }

  // Give PNG the file handle
  png_init_io(png_ptr, fp);

  // Set PNG Header
  png_set_IHDR(png_ptr, info_ptr, 
      width,                          // Width
      height,                         // Height
      8,                              // Bit depth
      PNG_COLOR_TYPE_RGB_ALPHA,       // Color type
      PNG_INTERLACE_NONE,             // Interlacing
      PNG_COMPRESSION_TYPE_DEFAULT,   // Compression
      PNG_FILTER_TYPE_DEFAULT);       // Filter method


  // Allocate a pointer to an array of png_byte pointers
  png_bytep * row_pointers = png_malloc(png_ptr, height * png_sizeof(png_bytep));

  // OpenGL stores pixel data in row major format, from bottom of image to top.
  // PNG's are stored row major, top to bottom.
  int i;
  for (i = 0; i < height; ++i)
    row_pointers[i] = &((png_bytep) buffer)[(height - i - 1) * width * 4];

  // Set the rows
  png_set_rows(png_ptr, info_ptr, row_pointers);
  // Write the png
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);


  free(row_pointers);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);

  return 0;
} // writepng()
