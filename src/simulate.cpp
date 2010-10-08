#include <iostream>
#include <fstream>
#include <string>

#include "slotted_discs.h"
#include "getopt.h"

void processOptions(string & name, SlottedDiscs * discs,
                    int argc, char ** argv);

int main(int argc, char ** argv)
{
  SlottedDiscs * discs = new SlottedDiscs();
  string outputFilename;
  
  processOptions(outputFilename, discs, argc, argv);
  ofstream OutputFile(outputFilename.c_str(), ios::binary);
  cout << "Slotted Discs C++ Simulation" << endl;
  cout << "Parameters" << endl;
  discs->printParameters();
  cout << "Initial Conditions" << endl;
  discs->printState();

  // Write initial condition record data
  OutputFile << discs;

  bool tensile[2] = {false, false};
  double tj, state[6] = {discs->q1, discs->q2, discs->q3,
                         discs->q4, discs->q5, discs->w};
  for (int j = 1; j < discs->fps*discs->tf + 1; ++j) {
    tj = ((double) j) / ((double) discs->fps);
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
  cout << "Simulation completed.  Data written to " << 
       outputFilename << endl;
  return 0;
} // main

void processOptions(string & name, SlottedDiscs * discs, int argc, char ** argv)
{
  // Default parameters
  DiscParams * p = new DiscParams;
  double ma, mb, ra, rb, l, g, alpha;
  ma = mb = 2.0;
  ra = rb = 0.1;
  l = sqrt(2.0)*ra;
  g = 9.81;
  alpha = M_PI/2.0;

  // Default initial states
  double state[6] = {discs->q1, discs->q2, discs->q3, discs->q4, discs->q5,
                     discs->w};

  // Default filename
  name = "./simulation.dat";

  int c, opt_index;
  struct option long_options[] = {
     {"help", no_argument, 0, '?'},
     {"ma",  required_argument, 0, 'm'},
     {"mb",  required_argument, 0, 'n'},
     {"ra",  required_argument, 0, 'o'},
     {"rb",  required_argument, 0, 'p'},
     {"alpha",  required_argument, 0, 'a'},
     {"offset",   required_argument, 0, 'l'},
     {"gravity", required_argument, 0, 'g'},
     {"heading",  required_argument, 0, 'h'},
     {"roll",  required_argument, 0, 'r'},
     {"spin",  required_argument, 0, 's'},
     {"x",  required_argument, 0, 'x'},
     {"y",  required_argument, 0, 'y'},
     {"omega",  required_argument, 0, 'w'},
     {"time",  required_argument, 0, 't'},
     {"fps",  required_argument, 0, 'f'},
     {"outputfile",  required_argument, 0, 'q'},
     {0, 0, 0, 0} };
  while (1) {
    opt_index = 0;
    c = getopt_long(argc, argv, "?m:n:o:p:a:l:g:h:r:s:x:y:w:t:f:q:", long_options, &opt_index);

  if (c == -1)
    break;

  switch (c) {
    case '?':
      printf(
"usage: %s [OPTION]\n\n"
"  -?, --help                  Display this help and exit\n"
"  --ma=val                    Mass of disc A\n"
"  --mb=val                    Mass of disc B\n"
"  --ra=val                    Radius of disc A\n"
"  --rb=val                    Radius of disc B\n"
"  --alpha=val                 Angle between disc planes\n"
"  --offset=val                Offset between disc centers\n"
"  --gravity=val               Gravity\n"
"  --heading=val               Initial heading of Disc A\n"
"  --roll=val                  Initial roll of Disc A\n"
"  --spin=val                  Initial spin of Disc A\n"
"  --x=val                     Initial x coordinate of Disc A contact\n"
"  --y=val                     Initial y coordinate of Disc A contact\n"
"  --omega=val                 Initial angular velocity about contact line\n"
"  --time=val                  Simulation time\n"
"  --fps=val                   Output data rate\n"
"  --outputfile=val            Output filename\n\n"
"Example of how to specify ma=1.0, mb=2.0, ra=0.1, intial angular velocity\n"
"of w=0.1:\n\n"
"$ %s --ma=1.0 --mb=2.0 --ra=0.1 --w=0.1\n\n", argv[0], argv[0]);
      exit(0);
 
    case 'm': ma = atof(optarg); break;
    case 'n': mb = atof(optarg); break;
    case 'o': ra = atof(optarg); break;
    case 'p': rb = atof(optarg); break;
    case 'a': alpha = atof(optarg); break;
    case 'l': l = atof(optarg); break;
    case 'g': g = atof(optarg); break;
    case 'h': state[0] = atof(optarg); break;
    case 'r': state[1] = atof(optarg); break;
    case 's': state[2] = atof(optarg); break;
    case 'x': state[3] = atof(optarg); break;
    case 'y': state[4] = atof(optarg); break;
    case 'w': state[5] = atof(optarg); break;
    case 't': discs->tf = atof(optarg); break;
    case 'f': discs->fps = atof(optarg); break;
    case 'q': name = optarg; break;
    default: abort();
    } // switch(c)
  } // while

  setParams(p, ma, mb, ra, rb, l, alpha, g);
  discs->setParameters(p);
  discs->setState(state);
  discs->evalConstants();
  discs->eoms();
  discs->computeOutputs();
} // processOptions()
