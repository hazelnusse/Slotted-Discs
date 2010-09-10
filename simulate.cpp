#include <iostream>
#include <fstream>

#include "twindiscs.h"

// Function forward declaration
void setParams(DiscParams * p, double ma, double mb, double ra, double rb, double l, double alpha, double g);

int main(int argc, char ** argv)
{
  SlottedDiscs * discs = new SlottedDiscs();
  ofstream OutputFile("./simulatecpp.dat", ios::binary);
  DiscParams * p = new DiscParams;
  
  double ma, mb, ra, rb, l, alpha, g;
  ma = mb = 2.0;
  ra = rb = 0.1;
  l = 0.15;//  sqrt(2.0)*ra;
  alpha = M_PI/2.0;
  g = 9.81;
  setParams(p, ma, mb, ra, rb, l, alpha, g);
  discs->setParameters(p);

  // Numerical integration loop
  int fps = 100;
  double tj, state[6] = {0.0, M_PI/4.0, M_PI/2.0, 0.5, 0.5, 1.0};
  // Set the speed if that is all that is specified
  if (argc == 2)
    state[5] = atof(argv[1]);
  // Set the speed and the simulation time if both are specified
  if (argc == 3) {
    state[5] = atof(argv[1]);
    discs->tf = atof(argv[2]);
  }

  discs->setState(state);
  discs->eoms();
  discs->computeOutputs();
  cout << "Parameters" << endl;
  discs->printParameters();
  cout << "Initial Conditions" << endl;
  discs->printState();

  // Write initial condition record data
  OutputFile << discs;

  for (int j = 1; j < fps*discs->tf + 1; ++j) {
    tj = ((double) j) / ((double) fps);
    while (discs->t < tj)
      gsl_odeiv_evolve_apply(discs->e, discs->c, discs->s,
                             &(discs->sys), &(discs->t), tj,
                             &(discs->h), state);
    discs->computeOutputs();
    OutputFile << discs;
  } // for i

  OutputFile.close();
  delete discs;
  delete p;
  return 0;
} // main

void setParams(DiscParams * p, double ma, double mb, double ra, double rb, double l, double alpha, double g) 
{
  double Ia = ma*ra*ra/4.0;
  double Ja = ma*ra*ra/2.0;
  double Ib = mb*rb*rb/4.0;
  double Jb = mb*rb*rb/2.0;
  p->m = ma + mb;
  p->ra = ra;
  p->rb = rb;
  p->l = l;
  p->g = 9.81;
  p->alpha = alpha;
  p->k = l*mb/(ma+mb);
  p->Ixx = Ia + Ib*pow(cos(alpha),2) + Jb*pow(sin(alpha),2) + mb*pow(l,2)*(ma*mb/
  pow((ma+mb),2)+pow((1-mb/(ma+mb)),2));
  p->Iyy = Ja + Jb + pow(sin(alpha),2)*(Ib-Jb) + mb*pow(l,2)*(ma*mb/pow((ma+mb),
  2)+pow((1-mb/(ma+mb)),2));
  p->Izz = Ia + Ib + mb*pow(l,2)*pow((1-mb/(ma+mb)),2) + mb*pow(l,2)*(-1+mb/(ma+
  mb))*(1-mb/(ma+mb));
  p->Ixy = sin(alpha)*cos(alpha)*(Ib-Jb);
} // setParams
