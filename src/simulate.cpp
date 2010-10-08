#include <iostream>
#include <fstream>

#include "slotted_discs.h"

int main(int argc, char ** argv)
{
  SlottedDiscs * discs = new SlottedDiscs();
  ofstream OutputFile("./simulation.dat", ios::binary);
  DiscParams * p = new DiscParams;
  
  double ma, mb, ra, rb, l, alpha, g;
  ma = mb = 2.0;
  ra = rb = 0.1;
  l = .1;//sqrt(2.0)*ra;
  alpha = M_PI/2.0;
  g = 9.81;
  setParams(p, ma, mb, ra, rb, l, alpha, g);
  discs->setParameters(p);

  // Numerical integration loop
  int fps = 1000;
  double tj, state[6] = {0.0, M_PI/4.0, M_PI/2.0, 0.0, 0.0, 1.0};
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

  bool tensile[2] = {false, false};
  for (int j = 1; j < fps*discs->tf + 1; ++j) {
    tj = ((double) j) / ((double) fps);
    while (discs->t < tj)
      gsl_odeiv_evolve_apply(discs->e, discs->c, discs->s,
                             &(discs->sys), &(discs->t), tj,
                             &(discs->h), state);
    discs->computeOutputs();
    if (discs->tensile[0] != tensile[0]) { // Disc A tensile force switched sign
      cerr << "Disc A tensile force sign change:  t = "
           << discs->t << ", faz = " << discs->faz << endl;
      tensile[0] = discs->tensile[0];
    }
    if (discs->tensile[1] != tensile[1]) { // Disc B tensile force switched sign
      cerr << "Disc B tensile force sign change:  t = "
           << discs->t << ", fbz = " << discs->fbz << endl;
      tensile[1] = discs->tensile[1];
    }
    OutputFile << discs;
  } // for j

  OutputFile.close();
  delete discs;
  delete p;
  return 0;
} // main
