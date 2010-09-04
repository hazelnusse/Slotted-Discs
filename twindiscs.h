#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

using namespace std;

extern "C" {
  inline int eomwrapper(double t, const double x[6], double f[6], void * params);
}

typedef struct {
  double ma, mb, ra, rb, l, g, alpha, Ia, Ib, Ja, Jb;
} DiscParams;

class SlottedDiscs {
  private:

  public:
    double t, tf, h, f[6], state[6];
    double q1, q2, q3, q4, q5, w1, w2, w3, q1p, q2p, q3p, q4p, q5p, w1p;
    double ma, mb, ra, rb, ke, pe;
    double Ia, Ib, Ja, Jb, g, alpha;
    double m, l, k, Ixx, Iyy, Izz, Iyz;
    double z[323];
    double T_da[16], T_db[16], T_ca[16], T_cb[16], T_so[16];
    double con[3], A[36], dfdt[6], p[3], H[3];

    // Camera settings
    double theta, phi, spin;
    double camx, camy, camz, dcam;
    
    // Numerical integrator variables
    const gsl_odeiv_step_type * T;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c; gsl_odeiv_evolve * e;
    gsl_odeiv_system sys;

    // Member functions
    SlottedDiscs(double * state = NULL, DiscParams * p = NULL);
    ~SlottedDiscs();

    // Accessors
    void printConstraints(void) const;
    void printEnergy(void) const;
    void printParameters() const;
    void printState(void) const;

    // Mutators
    int computeOutputs(void);
    int eoms(void);
    int setIntegrator(double t_i = 0.0,
                      double h_i = 0.001,
                      double t_f = 5.0,
                      const gsl_odeiv_step_type * int_T = gsl_odeiv_step_rk8pd);
    int setParams(DiscParams * p);
    int setState(double * x);
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);


}; 
