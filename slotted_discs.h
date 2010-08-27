#ifndef SLOTTED_DISCS_H
#define SLOTTED_DISCS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>

typedef struct {
  double t, tf, h, x[6], f[6], w1, w2;
  double Ix, Iy, Iz, g, alpha;
  double ma, mb, m, ra, rb, k, l, ke, pe;
  double z[222];
  double T_da[16], T_db[16];
  double no_cb[3];
  double A[4];
  const gsl_odeiv_step_type * T;
  gsl_odeiv_step * s;
  gsl_odeiv_control * c;
  gsl_odeiv_evolve * e;
  gsl_odeiv_system sys;
} sd_t;

void sdInit(sd_t * p);

int sdF(double t, const double *x, double *f, void *params);

void sdSetInertia(sd_t * p);

void sdOutputs(sd_t * p);

void sdPrint(sd_t * P);


#endif
