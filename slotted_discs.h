#ifndef SLOTTED_DISCS_H
#define SLOTTED_DISCS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>

typedef struct {
  double t, tf, h, x[5], f[5], w1, w2;
  double Ix, Iy, Iz, g, alpha;
  double ma, mb, m, ra, rb, k, l, ke, pe;
  double z[204];
  double T_da[16], T_db[16];
  double no_cb[3];
  double A[4];
  double q3, q3p;
  const gsl_odeiv_step_type * T;

  // Numerical integrator variables
  gsl_odeiv_step * s;
  gsl_odeiv_control * c;
  gsl_odeiv_evolve * e;
  gsl_odeiv_system sys;

  // Newton-Raphson variables
  gsl_root_fdfsolver * fdf_s;
  const gsl_root_fdfsolver_type * fdf_t;
  gsl_function_fdf fdf;
  int status;
} sd_t;

double sdcalcSpin(sd_t * p);

int sdF(double t, const double *x, double *f, void *params);

double sdConstraint_f(double q3, void * p);

double sdConstraint_df(double q3, void * p);

void sdConstraint_fdf(double q3, void * p, double * f, double * df);

void sdFree(sd_t *p);

void sdInit(sd_t * p);

void sdSetInertia(sd_t * p);

void sdOutputs(sd_t * p);

void sdPrint(sd_t * P);

#endif
