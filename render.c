#include <stdio.h>
#include <gsl/gsl_odeiv.h>
#include <GL/osmesa.h>
#include <GL/glu.h>


#include "writepng.h"
#include <png.h>
#include "slotted_discs.h"

// Forward declarations
static void render_image(sd_t * p, int Width, int Height);

int main(int argc, char ** argv)
{
  // OSMesa variables
  OSMesaContext ctx;
  void *buffer;
  char filename[255], *filename_root;
  int Width = 1920, Height = 1080;

  // Simulation variables
  sd_t * p = (sd_t *) malloc(sizeof(sd_t));
  sdInit(p);
  int i, fps = 60;
  double tj;

  // Checking for proper usage
  if (argc < 3) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s w3 filename_root [width height]\n\n", argv[0]);
    fprintf(stderr, "     w3 : initial angular rate about axis connecting discs.\n");
    fprintf(stderr, "     filename_root : Base file name for png files to render frames to.\n");
    return 0;
  }

  p->x[5] = atof(argv[1]);
  filename_root = argv[3];
  if (argc == 6) {
    p->tf = atof(argv[2]);
    Width = atoi(argv[4]);
    Height = atoi(argv[5]);
  }

  /* Create an RGBA-mode context */
  /* specify Z, stencil, accum bit sizes */
  ctx = OSMesaCreateContextExt(
                OSMESA_RGBA,  // RGBA Mode
                16,           // Z buffer bit depth
                0,            // stencil buffer bit depth
                0,            // accumulator buffer bit depth
                NULL);
  if (!ctx) {
    printf("OSMesaCreateContext failed!\n");
    return 0;
  }

  /* Allocate the image buffer */
  buffer = malloc(Width * Height * 4 * sizeof(GLubyte));
  if (!buffer) {
    printf("Alloc image buffer failed!\n");
    return 0;
  }

  /* Bind the buffer to the context and make it current */
  if (!OSMesaMakeCurrent(ctx, buffer, GL_UNSIGNED_BYTE, Width, Height)) {
    printf("OSMesaMakeCurrent failed!\n");
    return 0;
  }
  
  sdF(p->t, p->x, p->f, p);
  sdOutputs(p);
  render_image(p, Width, Height);
  sprintf(filename, "%s000.png", filename_root);
  writepng(buffer, filename, Width, Height);
  for (i = 1; i < fps*p->tf + 1; ++i) {
    tj = ((double) i) / ((double) fps);
    while (p->t < tj) {
      gsl_odeiv_evolve_apply(p->e, p->c, p->s,
          &(p->sys), &(p->t), tj, &(p->h), p->x);
    }
    sdOutputs(p);
    render_image(p, Width, Height);
    sprintf(filename, "%s%03d.png", filename_root, i);
    writepng(buffer, filename, Width, Height);
  } // for i

  sdFree(p);
  // free the image buffer
  free(buffer);

  // destroy the context
  OSMesaDestroyContext(ctx);

  return 0;
} // main()


static void Disc(double r)
{
  GLUquadric *q = gluNewQuadric();
  gluQuadricDrawStyle(q, GLU_FILL);
  gluQuadricNormals(q, GLU_SMOOTH);
  // Drawing routines
  glPushMatrix();
    // Rotate about x so cylinder height is along the previous y axis
    glRotatef(90.0, 1.0, 0.0, 0.0);
    gluCylinder(q, r, 0.0, .1*r, 100, 50);
    glPushMatrix();
      glRotatef(180.0, 1.0, 0.0, 0.0);
      gluCylinder(q, r, 0.0, .1*r, 100, 50);
    glPopMatrix();
  glPopMatrix();
  gluDeleteQuadric(q);
} // Disc()
/*
static void ReferenceFrame(void)
{
  GLUquadric *q = gluNewQuadric();
  gluQuadricDrawStyle(q, GLU_FILL);
  gluQuadricNormals(q, GLU_SMOOTH);
  gluCylinder(q, .01, .01, .05, 10, 10);
  gluDeleteQuadric(q);
}  */

static void render_image(sd_t * p, int Width, int Height)
{
   GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
   GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
   GLfloat red_mat[]   = { 1.0, 0.2, 0.2, 1.0 };
   GLfloat green_mat[] = { 0.2, 1.0, 0.2, 1.0 };
   // GLfloat blue_mat[]  = { 0.2, 0.2, 1.0, 1.0 };


   glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
   glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   
   glClearColor(0.0, 0.0, 0.0, 1.0);
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
   
   glViewport(0, 0, (GLsizei) Width, (GLsizei) Height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(45.0, (GLfloat) Width/(GLfloat) Height, .5, 5.0);

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   //glRotatef(90.0, 1.0, 0.0, 0.0);
   //glRotatef(-90.0, 0.0, 0.0, 1.0);
   gluLookAt(p->x[3] + 0.5, p->x[4], -1.0,   // camera position 
             0.0, 0.0, 0.0,   // point camera at this position 
             0.0, 0.0, -1.0);  // define up of the camera 

   glPushMatrix();

   glPushMatrix();
     glMultMatrixd(p->T_da);
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red_mat);
     Disc(p->ra);
   glPopMatrix();
   
   glPushMatrix();
     glMultMatrixd(p->T_db);
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, green_mat);
     Disc(p->rb);
   glPopMatrix();


   glPopMatrix();

   glFinish();
} // render_image()
