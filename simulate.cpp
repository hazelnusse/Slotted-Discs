#include <iostream>
#include "twindiscs.h"

int main(int argc, char ** argv)
{
  SlottedDiscs * discs = new SlottedDiscs();
  discs->eoms();
  discs->computeOutputs();
  cout << "State" << endl;
  discs->printState();
  cout << "Constraints" << endl;
  discs->printConstraints();
  cout << "Energy" << endl;
  discs->printEnergy();
  cout << "Parameters" << endl;
  discs->printParameters();
  delete discs;
  return 0;
} // main
