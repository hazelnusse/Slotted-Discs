#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <Eigen/Eigen>

#define Z_MAX 386

using namespace std;

extern "C" {
  inline int eomwrapper(double t, const double x[6], double f[6], void * params);
}

typedef struct {
  double m, ra, rb, l, k, g, alpha, Ixx, Iyy, Izz, Ixy;
} DiscParams;

typedef Eigen::Matrix<double, 3, 3> Matrix33d;

// Given masses, radii, offset distance and angles of two discs, populate a
// DiscParams struct with the minimal model parameters.
void setParams(DiscParams * p, double ma, double mb, double ra, double rb, double l, double alpha, double g);

class SlottedDiscs {
  private:

  public:
    double ra, rb, l;
    double Ia, Ib, Ja, Jb, g, alpha;
    double m, k, Ixx, Iyy, Izz, Ixy;

    double t, tf, h;
    double q1, q2, q3, q4, q5, w, w1, w2, w3, q1p, q2p, q3p, q4p, q5p, wp;
    double ke, pe, te, equilibria;
    double z[Z_MAX];
    double no_cb[3], no_so[3], a_so[3], H[3], p[3], df[36];
    double T_da[16], T_db[16], T_so[16], T_ca[16], T_cb[16], T_dagl[16], T_dbgl[16];
    double fx, fay, faz, fby, fbz;

    // Variables for holonomic constraint
    double hc_f, hc_df[2];

    // Boolean for deciding whether to integrate dq2/dt or dq3/dt
    // integrateq2dot == True --> q3 is dependent coordinate
    // integrateq2dot == False --> q2 is dependent coordinate
    bool integrateq2dot;

    // Camera variables
    double theta, phi, d, ctx, cty, ctz;

    // Flag array to detect when forces are tensile
    bool tensile[2];

    // Eigenvalue variables
    Matrix33d A;

    // Numerical integrator variables
    const gsl_odeiv_step_type * T;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c; gsl_odeiv_evolve * e;
    gsl_odeiv_system sys;
    int fps;

    // Member functions
    SlottedDiscs();
    ~SlottedDiscs();

    // Accessors
    void printEnergy(void) const;
    void printParameters() const;
    void printState(void) const;
    void writeRecord_dt(void) const;
    double hc(double q2, double q3) const;
    double hc_dq2(double q2, double q3) const;
    double hc_dq3(double q2, double q3) const;

    // Mutators
    void setState(const double state[6]);
    void setParameters(DiscParams * p);
    void evalConstants(void);
    void eoms(void);
    void computeOutputs(void);

    // Wrapper functions to interface with GSL required calling conventions
    friend int eomwrapper(double t, const double x[6], double f[6], void * params);
    friend ostream &operator<<(ostream &file, const SlottedDiscs *discs);
};
