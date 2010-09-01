#include <stdio.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include "slotted_discs.h"

int main(int argc, char ** argv)
{
  sd_t * p = (sd_t *) malloc(sizeof(sd_t));
  int i, fps = 100;
  double tj;
  FILE * fp = fopen("./simulate.dat", "wb");
  

  sdInit(p);
  sdWriteRecord(p, fp);
  p->tf = 20.0;
  for (i = 1; i < fps*p->tf + 1; ++i) {
    tj = ((double) i) / ((double) fps);
    while (p->t < tj)
      gsl_odeiv_evolve_apply(p->e, p->c, p->s,
          &(p->sys), &(p->t), tj, &(p->h), p->x);
    sdOutputs(p);
    sdWriteRecord(p, fp);
  } // for i

  sdFree(p);
  fclose(fp);
  return 0;
} // main
