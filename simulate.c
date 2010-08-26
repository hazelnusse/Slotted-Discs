#include <stdio.h>
#include <gsl/gsl_odeiv.h>
#include "slotted_discs.h"

int main(int argc, char ** argv)
{
  sd_t * discs = (sd_t *) malloc(sizeof(sd_t));
  sdInit(discs);
  int i, tj;

  for (i = 1; i < 60*discs->tf + 1; ++i) {
    tj = i / 60.0;
    while (discs->t < tj)
      gsl_odeiv_evolve_apply(discs->e, discs->c, discs->s, &(discs->sys),
          &(discs->t), tj, &(discs->h), discs->x);
  } // for i


  return 0;
} // main
