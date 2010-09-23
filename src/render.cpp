#include <iostream>
// #include <png.h>
#include <gsl/gsl_odeiv.h>
#include <GL/osmesa.h>
#include <GL/glu.h>


#include "writepng.h"
#include "slotted_discs.h"

// Forward declarations
static void render_image(SlottedDiscs * p, int Width, int Height);

int main(int argc, char ** argv)
{
  // OSMesa variables
  OSMesaContext ctx;
  void *buffer;
  char filename[255], *filename_root;
  int Width = 1920, Height = 1080;

  // Simulation variables
  SlottedDiscs * discs = new SlottedDiscs();
  DiscParams * p = new DiscParams;
  int fps = 60;
  double tj;
  double state[6] = {M_PI, M_PI/4.0, M_PI/2.0, 0.0, 0.0, 1.0};
  
  double ma, mb, ra, rb, l, alpha, g;
  ma = mb = 2.0;
  ra = rb = 0.1;
  l = .1; //sqrt(2.0)*ra;
  alpha = M_PI/2.0;
  g = 9.81;
  setParams(p, ma, mb, ra, rb, l, alpha, g);
  discs->setParameters(p);

  // Checking for proper usage
  if (argc < 3) {
    cerr << "Usage:" << endl;
    cerr << "  %s" << argv[0] << " w tf filename_root [width height]" << endl << endl;
    cerr << "     w : initial angular rate about contact line." << endl;
    cerr << "     tf : simulation time." << endl;
    cerr << "     filename_root : Base file name for png files to render frames to." << endl;
    return 0;
  }

  state[5] = atof(argv[1]);
  discs->tf = atof(argv[2]);
  discs->d = atof(argv[3]);
  filename_root = argv[4];
  if (argc == 7) {
    Width = atoi(argv[5]);
    Height = atoi(argv[6]);
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
  
  discs->setState(state);
  discs->eoms();
  discs->computeOutputs();
  render_image(discs, Width, Height);
  sprintf(filename, "%s0000.png", filename_root);
  writepng(buffer, filename, Width, Height);
  for (int j = 1; j < fps*discs->tf + 1; ++j) {
    tj = ((double) j) / ((double) fps);
    while (discs->t < tj)
      gsl_odeiv_evolve_apply(discs->e, discs->c, discs->s,
                             &(discs->sys), &(discs->t), tj,
                             &(discs->h), state);
    discs->computeOutputs();
    render_image(discs, Width, Height);
    sprintf(filename, "%s%04d.png", filename_root, j);
    writepng(buffer, filename, Width, Height);
  } // for j

  // free the image buffer
  delete discs;
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
    gluCylinder(q, r, 0.0, 0.2*r, 100, 50);
    glRotatef(180.0, 1.0, 0.0, 0.0);
    gluCylinder(q, r, 0.0, 0.2*r, 100, 50);
  glPopMatrix();
  gluDeleteQuadric(q);
} // Disc()
/*static void ReferenceFrame(void)
{
  GLUquadric *q = gluNewQuadric();
  gluQuadricDrawStyle(q, GLU_FILL);
  gluQuadricNormals(q, GLU_SMOOTH);
  gluCylinder(q, .01, .01, .05, 10, 10);
  gluDeleteQuadric(q);
}*/

static void render_image(SlottedDiscs * p, int Width, int Height)
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

   glPushMatrix();
     glMultMatrixd(p->T_dagl);
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red_mat);
     Disc(p->ra);
   glPopMatrix();
   
   glPushMatrix();
     glMultMatrixd(p->T_dbgl);
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, green_mat);
     Disc(p->rb);
   glPopMatrix();

   glFinish();
} // render_image()
