#include <stdio.h>
#include <gsl/gsl_odeiv.h>
#include "slotted_discs.h"

int main(int argc, char ** argv)
{
  sd_t * discs = (sd_t *) malloc(sizeof(sd_t));
  sdInit(discs);
  sdF(0.0, discs->x, discs->f, discs);
  sdOutputs(discs);
  discs->x[5] = 1.0;
  int i, fps = 10;
  double tj;

  for (i = 0; i < 222; ++i)
    printf("discs->z[%d] = %0.5g\n", i, discs->z[i]);

  printf("%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s\n",
         "t", "yaw", "roll", "spin", "x", "y", "u", "nocbx", "nocby", "nocbz");
  sdPrint(discs);
  for (i = 1; i < fps*discs->tf + 1; ++i) {
    tj = ((double) i) / ((double) fps);
    while (discs->t < tj) {
      gsl_odeiv_evolve_apply(discs->e, discs->c, discs->s,
          &(discs->sys), &(discs->t), tj, &(discs->h), discs->x);
    
      sdPrint(discs);
    }
  } // for i


  gsl_odeiv_evolve_free(discs->e);
  gsl_odeiv_control_free(discs->c);
  gsl_odeiv_step_free(discs->s);
  return 0;
} // main
