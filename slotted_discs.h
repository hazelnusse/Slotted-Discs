#ifndef SLOTTED_DISCS_H
#define SLOTTED_DISCS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>
#include <math.h>

typedef struct {
  double t, tf, h, x[6], f[6], w1, w2;
  double Ia, Ib, Ja, Jb, g, alpha;
  double ma, mb, ra, rb, l, ke, pe;
  double z[383];
  double T_da[16], T_db[16];
  double no_cb[3];
  double A[36];
  double dfdt[6];

  // Numerical integrator variables
  const gsl_odeiv_step_type * T;
  gsl_odeiv_step * s;
  gsl_odeiv_control * c; gsl_odeiv_evolve * e;
  gsl_odeiv_system sys;

  // Newton-Raphson variables
  gsl_root_fdfsolver * fdf_s;
  const gsl_root_fdfsolver_type * fdf_t;
  gsl_function_fdf fdf;
  int status;
} sd_t;

double sdCalcSpin(sd_t * p);

double sdConstraint_f(double q3, void * p);

double sdConstraint_df(double q3, void * p);

void sdConstraint_fdf(double q3, void * p, double * f, double * df);

int sdF(double t, const double * x, double * f, void * params);

int sdJ(double t, const double * x, double * dfdx, double *dfdt, void * params);

void sdFree(sd_t *p);

void sdInit(sd_t * p);

void sdInitOdeSolver(sd_t * p);

void sdInitParameters(sd_t * p);

void sdInitRootFinder(sd_t * p);

void sdInitState(sd_t * p);

void sdOutputs(sd_t * p);

void sdPrint(sd_t * P);

void sdPrintf(sd_t * P);

void sdWriteRecord(sd_t * p, FILE * fp);

#endif
