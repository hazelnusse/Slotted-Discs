#include <stdio.h>
#include <gsl/gsl_odeiv.h>
#include "slotted_discs.h"

int main(int argc, char ** argv)
{
  sd_t * p = (sd_t *) malloc(sizeof(sd_t));
  int i, fps = 100;
  double tj, w3;
  FILE * fp = fopen("./simulate.dat", "wb");
  
  if (argc == 2)
    w3 = atof(*(argv + 1));


  sdInit(p);
  p->x[5] = w3;
  sdF(p->t, p->x, p->f, p);
  sdOutputs(p);
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
