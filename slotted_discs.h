#ifndef SLOTTED_DISCS_H
#define SLOTTED_DISCS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>

typedef struct {
  double t, x[6], f[6];
  double Ix, Iy, Iz, g, alpha;
  double m, ra, rb, k, l, ke, pe;
  double z[222];
  double T_da[16], T_db[16];
  double no_cb[3];
  double A[4];
} sd_t;

void sdInit(sd_t * p);

int sdF(double t, const double *x, double *f, void *params);

int sdOutputs(sd_t * p);

#endif
