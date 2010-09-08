#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

using namespace std;

extern "C" {
  int eomwrapper(double t, const double x[6], double f[6], void * params);
}

typedef struct {
  double ma, mb, ra, rb, l, g, alpha;
} DiscParams;

class SlottedDiscs {
  private:

  public:
    double t, tf, h, f[6];
    double q1, q2, q3, q4, q5, w1, w2, w3, q1p, q2p, q3p, q4p, q5p, w3p;
    double ma, mb, ra, rb, l;
    double ke, pe, te;
    double Ia, Ib, Ja, Jb, g, alpha;
    // double m, l, k, Ixx, Iyy, Izz, Iyz;
    double z[171];
    double no_cb[3];

    // Numerical integrator variables
    const gsl_odeiv_step_type * T;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c; gsl_odeiv_evolve * e;
    gsl_odeiv_system sys;

    // Member functions
    SlottedDiscs();
    ~SlottedDiscs();

    // Accessors
    void printConstraints(void) const;
    void printEnergy(void) const;
    void printParameters() const;
    void printState(void) const;
    friend ostream &operator<<(ostream &file, const SlottedDiscs *discs);

    // Mutators
    void setState(double state[6]);
    void setParameters(DiscParams * p);
    void computeOutputs(void);
    int eoms(void);
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);
}; 
