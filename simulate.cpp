#include <iostream>
#include <fstream>

#include "twindiscs.h"

int main(int argc, char ** argv)
{
  SlottedDiscs * discs = new SlottedDiscs();
  ofstream OutputFile("./simulatecpp.dat", ios::binary);
  
  DiscParams p;
  p.ma = p.mb = 2.0;
  p.ra = p.rb = 0.1;
  p.l = 0.15;//  sqrt(2.0)*p.ra;
  p.g = 9.81;
  p.alpha = M_PI/2.0;
  discs->setParameters(&p);

  // Numerical integration loop
  int fps = 100;
  double tj, state[6] = {0.0, M_PI/4.0, M_PI/2.0, 0.5, 0.5, 1.0};
  if (argc == 2)
    state[5] = atof(argv[1]);
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
  return 0;
} // main
