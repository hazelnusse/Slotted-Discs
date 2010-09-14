#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#define Z_MAX 134

using namespace std;

extern "C" {
  int eomwrapper(double t, const double x[6], double f[6], void * params);
}

typedef struct {
  double m, ra, rb, l, k, g, alpha, Ixx, Iyy, Izz, Ixy;
} DiscParams;

class SlottedDiscs {
  private:

  public:
    double ra, rb, l;
    double Ia, Ib, Ja, Jb, g, alpha;
    double m, k, Ixx, Iyy, Izz, Ixy;
    
    double t, tf, h;
    double q1, q2, q3, q4, q5, w1, w2, w3, q1p, q2p, q3p, q4p, q5p, w3p;
    double ke, pe, te;
    double z[Z_MAX];
    double no_cb[2], con[3];

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
    void writeRecord_dt(void) const;
    friend ostream &operator<<(ostream &file, const SlottedDiscs *discs);

    // Mutators
    void setState(const double state[6]);
    void setParameters(DiscParams * p);
    void evalConstants(void);
    void eoms(void);
    void computeOutputs(void);
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);
}; 
