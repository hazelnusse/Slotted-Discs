#include "slotted_discs.h"

int eomwrapper(double t, const double x[6], double f[6], void * params)
{
  SlottedDiscs * p = (SlottedDiscs *) params;
  // Assign the states of the SlottedDiscs object
  p->setState(x);
  // Evaluate the RHS of the ODE's representing the equations of motion
  p->eoms();
  // Assign the right hand sides of the ODE's to the array passed
  f[0] = p->q1p;
  f[1] = p->q2p;
  f[2] = p->q3p;
  f[3] = p->q4p;
  f[4] = p->q5p;
  f[5] = p->wp;
  // Return the status
  return GSL_SUCCESS;
}

ostream &operator<<(ostream &file, const SlottedDiscs * discs)
{
  file.write((char *) &(discs->t), sizeof discs->t);
  file.write((char *) &discs->q1, sizeof discs->q1);
  file.write((char *) &discs->q2, sizeof discs->q2);
  file.write((char *) &discs->q3, sizeof discs->q3);
  file.write((char *) &discs->q4, sizeof discs->q4);
  file.write((char *) &discs->q5, sizeof discs->q5);
  file.write((char *) &discs->w, sizeof discs->w);
  file.write((char *) &discs->w1, sizeof discs->w1);
  file.write((char *) &discs->w2, sizeof discs->w2);
  file.write((char *) &discs->w3, sizeof discs->w3);
  file.write((char *) &(discs->no_cb), sizeof discs->no_cb);
  file.write((char *) &(discs->ke), sizeof discs->ke);
  file.write((char *) &(discs->pe), sizeof discs->pe);
  file.write((char *) &(discs->te), sizeof discs->te);
  return file;
} // operator <<

void SlottedDiscs::writeRecord_dt(void) const
{
  ofstream fp("./record.py", ios::out);
  fp << "#!/usr/bin/env python" << endl;
  fp << "import numpy as np" << endl;
  fp << "record_dt = np.dtype([('t', np.float64), " <<
        "('q1', np.float64), " <<
        "('q2', np.float64), " <<
        "('q3', np.float64), " <<
        "('x', np.float64), " <<
        "('y', np.float64), " <<
        "('w', np.float64), " <<
        "('w1', np.float64), " <<
        "('w2', np.float64), " <<
        "('w3', np.float64), " <<
        "('cbx', np.float64), " <<
        "('cby', np.float64), " <<
        "('cbz', np.float64), " <<
        "('ke', np.float64), " <<
        "('pe', np.float64), " <<
        "('te', np.float64)]) " << endl;
  fp.close();
} // writeRecord_dt()

SlottedDiscs::SlottedDiscs()
{
  // Default parameters
  ra = rb = .1;
  l = sqrt(2.0)*ra;
  g = 9.81;
  alpha = M_PI/2.0;
  
  double ma, mb;
  ma = mb = 2.0;
  m = ma + mb;
  double Ia = ma*ra*ra/4.0;
  double Ja = ma*ra*ra/2.0;

  double Ib = mb*rb*rb/4.0;
  double Jb = mb*rb*rb/2.0;

  k = l*mb/(ma+mb);
  Ixx = Ia + Ib*pow(cos(alpha),2) + Jb*pow(sin(alpha),2) + mb*pow(l,2)*(ma*mb/
  pow((ma+mb),2)+pow((1-mb/(ma+mb)),2));
  Iyy = Ja + Jb + pow(sin(alpha),2)*(Ib-Jb) + mb*pow(l,2)*(ma*mb/pow((ma+mb),
  2)+pow((1-mb/(ma+mb)),2));
  Izz = Ia + Ib + mb*pow(l,2)*pow((1-mb/(ma+mb)),2) + mb*pow(l,2)*(-1+mb/(ma+
  mb))*(1-mb/(ma+mb));
  Ixy = sin(alpha)*cos(alpha)*(Ib-Jb);

  // Set state
  q1 = 0.0;
  q2 = M_PI / 4.0;
  q3 = M_PI/2.0;
  q4 = q5 = 0.0;
  w = -1.4;

  // Set integration settings
  t = 0.0;
  tf = 10.0;
  h = 0.001;
  T = gsl_odeiv_step_rk8pd;
  s = gsl_odeiv_step_alloc(T, 6);
  c = gsl_odeiv_control_y_new(1e-6, 1e-9);
  e = gsl_odeiv_evolve_alloc(6);
  sys.function = eomwrapper;
  sys.jacobian = NULL;
  sys.dimension = 6;
  sys.params = this;

  for (int i = 0; i < Z_MAX; ++i)
    z[i] = 0.0;

  // Constants
  evalConstants();
  eoms();
  computeOutputs();

  // Write the file which holds the numpy record data type so that plotting is
  // easy
  writeRecord_dt();
 
} // constructor

SlottedDiscs::~SlottedDiscs()
{
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
} // destructor

void SlottedDiscs::printState(void) const
{
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::scientific, ios::floatfield);
  cout.setf(ios::adjustfield, ios::right);
  cout.precision(9);
  cout.width(18);
  cout << q1;
  cout.width(18);
  cout << q2;
  cout.width(18);
  cout << q3;
  cout.width(18);
  cout << q4;
  cout.width(18);
  cout << q5;
  cout.width(18);
  cout << w3;
  cout << endl;
} // printState()

void SlottedDiscs::printEnergy(void) const
{
  cout.setf(ios::scientific, ios::floatfield);
  cout.setf(ios::adjustfield, ios::right);
  cout.precision(9);
  cout.width(18);
  cout << ke;
  cout.width(18);
  cout << pe;
  cout.width(18);
  cout << te;
  cout << endl;
} // printEnergy()

void SlottedDiscs::printParameters() const
{
  cout.setf(ios::scientific, ios::floatfield);
  cout.setf(ios::adjustfield, ios::right);
  cout.precision(9);
  cout << "ra = " << ra << endl << "rb = " << rb << endl;
  cout << "m = " << m << endl;
  cout << "Ixx = " << Ixx << endl << "Iyy = " << Iyy << endl;
  cout << "Izz = " << Izz << endl << "Ixy = " << Ixy << endl;
  cout << "alpha = " << alpha << endl << "l = " << l << endl;
  cout << "k = " << k << endl;
  cout << "g = " << g << endl;
} // printParameters()

void SlottedDiscs::evalConstants(void)
{
  z[97] = g*m;
  z[7] = cos(alpha);
  z[8] = sin(alpha);
} //evalConstants()

void SlottedDiscs::eoms(void)
{
  z[6] = sin(q3);
  z[9] = ra*z[6];
  z[4] = sin(q2);
  z[3] = cos(q2);
  z[13] = z[3]*z[6];
  z[20] = z[7]*z[4] + z[8]*z[13];
  z[21] = 1 - pow(z[20],2);
  z[22] = pow(z[21],0.5);
  z[24] = z[20]/z[22];
  z[23] = 1/z[22];
  z[25] = rb*(z[8]*z[24]-z[13]*z[23]);
  z[28] = z[9] + z[25];
  z[26] = rb*(z[7]*z[24]-z[4]*z[23]);
  z[5] = cos(q3);
  z[14] = z[3]*z[5];
  z[27] = rb*z[14]*z[23];
  z[10] = ra*z[5];
  z[29] = z[27] - l - z[10];
  z[30] = pow(z[26],2) + pow(z[28],2) + pow(z[29],2);
  z[31] = pow(z[30],0.5);
  z[32] = z[28]/z[31];
  w1 = w*z[32];
  z[33] = z[29]/z[31];
  w3 = w*z[33];
  q1p = -(sin(q3)*w1-cos(q3)*w3)/cos(q2);
  q2p = sin(q3)*w3 + cos(q3)*w1;
  z[11] = z[4]*z[6];
  z[18] = z[7]*z[3] - z[8]*z[11];
  z[36] = z[18]*z[20];
  z[38] = z[36]/pow(z[21],0.5);
  z[40] = (z[18]*z[22]+z[20]*z[38])/pow(z[22],2);
  z[42] = z[38]/pow(z[22],2);
  z[48] = rb*(z[7]*z[40]-z[3]*z[23]-z[4]*z[42]);
  z[44] = rb*(z[8]*z[40]+z[11]*z[23]-z[13]*z[42]);
  z[12] = z[4]*z[5];
  z[46] = rb*(z[12]*z[23]-z[14]*z[42]);
  z[51] = z[3]*z[26] + z[4]*z[48] + z[13]*z[44] + z[14]*z[46] - z[11]*(z[9]+
  z[25]) - z[12]*(l+z[10]-z[27]);
  z[35] = z[8]*z[14];
  z[37] = z[20]*z[35];
  z[39] = z[37]/pow(z[21],0.5);
  z[41] = (z[20]*z[39]+z[22]*z[35])/pow(z[22],2);
  z[43] = z[39]/pow(z[22],2);
  z[49] = rb*(z[7]*z[41]-z[4]*z[43]);
  z[45] = rb*(z[8]*z[41]-z[13]*z[43]-z[14]*z[23]);
  z[47] = rb*(z[13]*z[23]-z[14]*z[43]);
  z[50] = z[13]*(l-z[27]) - z[4]*z[49] - z[13]*z[45] - z[14]*z[25] - z[14]*
  z[47];
  z[52] = z[51]/z[50];
  q3p = z[52]*q2p;
  z[1] = cos(q1);
  z[2] = sin(q1);
  q4p = -(z[9]*(z[1]*z[6]+z[2]*z[12])+z[10]*(z[1]*z[5]-z[2]*z[11]))*q3p;
  q5p = -(z[10]*(z[1]*z[11]+z[2]*z[5])-z[9]*(z[1]*z[12]-z[2]*z[6]))*q3p;
  z[34] = z[26]/z[31];
  z[53] = z[34]*(k+z[10]);
  z[54] = z[9]*z[33] + z[32]*(k+z[10]);
  z[55] = z[9]*z[34];
  z[59] = 2*z[26]*z[48] + 2*z[28]*z[44] - 2*z[29]*z[46];
  z[61] = z[59]/pow(z[30],0.5);
  z[71] = (z[29]*z[61]+2*z[31]*z[46])/pow(z[31],2);
  z[73] = w*z[71];
  z[57] = z[10] + z[45];
  z[58] = z[9] - z[47];
  z[60] = 2*z[26]*z[49] + 2*z[28]*z[57] + 2*z[29]*z[58];
  z[62] = z[60]/pow(z[30],0.5);
  z[72] = (z[29]*z[62]-2*z[31]*z[58])/pow(z[31],2);
  z[74] = w*z[72];
  z[77] = -0.5*z[73]*q2p - 0.5*z[74]*q3p;
  z[63] = (z[28]*z[61]-2*z[31]*z[44])/pow(z[31],2);
  z[65] = w*z[63];
  z[64] = (z[28]*z[62]-2*z[31]*z[57])/pow(z[31],2);
  z[66] = w*z[64];
  z[75] = -0.5*z[65]*q2p - 0.5*z[66]*q3p;
  z[67] = (z[26]*z[61]-2*z[31]*z[48])/pow(z[31],2);
  z[69] = w*z[67];
  z[68] = (z[26]*z[62]-2*z[31]*z[49])/pow(z[31],2);
  z[70] = w*z[68];
  z[76] = 0.5*z[69]*q2p + 0.5*z[70]*q3p;
  z[81] = -z[9]*z[34] - 0.5*z[68]*(k+z[10]);
  z[83] = w*z[81];
  z[82] = z[67]*(k+z[10]);
  z[84] = w*z[82];
  z[93] = z[83]*q3p - pow(w,2)*z[33]*z[54] - pow(w,2)*z[34]*z[55] - 0.5*z[84]*
  q2p;
  z[85] = z[10]*z[33] - z[9]*z[32] - 0.5*z[9]*z[72] - 0.5*z[64]*(k+z[10]);
  z[87] = w*z[85];
  z[86] = -0.5*z[9]*z[71] - 0.5*z[63]*(k+z[10]);
  z[88] = w*z[86];
  z[94] = pow(w,2)*z[33]*z[53] + z[87]*q3p + z[88]*q2p - pow(w,2)*z[32]*z[55];
  z[89] = z[10]*z[34] - 0.5*z[9]*z[68];
  z[91] = w*z[89];
  z[90] = z[9]*z[67];
  z[92] = w*z[90];
  z[95] = pow(w,2)*z[32]*z[54] + pow(w,2)*z[34]*z[53] + z[91]*q3p - 0.5*z[92]*
  q2p;
  z[98] = z[97]*(z[13]*z[53]-z[4]*z[54]-z[14]*z[55]) + Izz*z[33]*z[77] + 
  z[32]*(Ixx*z[75]+Ixy*z[76]) + m*(z[53]*z[93]+z[54]*z[94]+z[55]*z[95]) - 
  z[34]*(Ixy*z[75]+Iyy*z[76]);
  z[96] = -Iyy*pow(z[34],2) - Izz*pow(z[33],2) - z[32]*(Ixx*z[32]-2*Ixy*z[34]) - 
  m*(pow(z[53],2)+pow(z[54],2)+pow(z[55],2));
  z[99] = z[98]/z[96];
  wp = z[99];
} // eoms()

void SlottedDiscs::computeOutputs(void)
{
  w2 = -w*z[34];
  ke = 0.5*pow(w,2)*(Iyy*pow(z[34],2)+Izz*pow(z[33],2)+z[32]*(Ixx*z[32]-2*Ixy*
  z[34])+m*(pow(z[53],2)+pow(z[54],2)+pow(z[55],2)));
  pe = g*m*(z[9]*z[13]+z[14]*(k+z[10]));
  te = ke + pe;
  z[100] = z[1]*z[5] - z[2]*z[11];
  z[101] = z[1]*z[11] + z[2]*z[5];
  z[102] = z[2]*z[3];
  z[103] = z[1]*z[3];
  z[104] = z[1]*z[6] + z[2]*z[12];
  z[105] = z[2]*z[6] - z[1]*z[12];

  no_cb[0] = q4 + z[26]*z[102] + z[100]*(z[9]+z[25]) - z[104]*(l+z[10]-z[27]);
  no_cb[1] = q5 + z[101]*(z[9]+z[25]) - z[26]*z[103] - z[105]*(l+z[10]-z[27]);
  no_cb[2] = -z[4]*z[26] - z[13]*(z[9]+z[25]) - z[14]*(l+z[10]-z[27]);
} // computeOutputs()

void SlottedDiscs::setState(const double state[6])
{
  q1 = state[0];
  q2 = state[1];
  q3 = state[2];
  q4 = state[3];
  q5 = state[4];
  w = state[5];
} // setState()

void SlottedDiscs::setParameters(DiscParams * p)
{
  m = p->m;
  ra = p->ra;
  rb = p->rb;
  k = p->k;
  l = p->l;
  g = p->g;
  alpha = p->alpha;
  Ixx = p->Ixx;
  Iyy = p->Iyy;
  Izz = p->Izz;
  Ixy = p->Ixy;

  evalConstants();
} // setParameters()
