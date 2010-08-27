#include <stdio.h>
#include <gsl/gsl_odeiv.h>
#include "slotted_discs.h"

int main(int argc, char ** argv)
{
  sd_t * discs = (sd_t *) malloc(sizeof(sd_t));
  sdInit(discs);
  discs->x[5] = 1.0;
  sdF(0.0, discs->x, discs->f, discs);
  sdOutputs(discs);
  int i, fps = 10;
  double tj;

  printf("%5s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n",
         "t", "yaw", "roll", "spin", "x", "y", "u", "nocbx",
         "nocby", "nocbz", "ke", "pe", "te");
  sdPrint(discs);
  for (i = 1; i < fps*discs->tf + 1; ++i) {
    tj = ((double) i) / ((double) fps);
    while (discs->t < tj)
      gsl_odeiv_evolve_apply(discs->e, discs->c, discs->s,
          &(discs->sys), &(discs->t), tj, &(discs->h), discs->x);
    sdOutputs(discs);
    sdPrint(discs);
  } // for i

  sdFree(discs);
  return 0;
} // main
