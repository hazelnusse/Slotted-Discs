#include "twindiscs.h"

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
  f[5] = p->w3p;
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
  file.write((char *) &discs->w1, sizeof discs->w1);
  file.write((char *) &discs->w2, sizeof discs->w2);
  file.write((char *) &discs->w3, sizeof discs->w3);
  file.write((char *) &(discs->con), sizeof discs->con);
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
        "('w1', np.float64), " <<
        "('w2', np.float64), " <<
        "('w3', np.float64), " <<
        "('nh1', np.float64), " <<
        "('nh2', np.float64), " <<
        "('hc', np.float64), " <<
        "('cbx', np.float64), " <<
        "('cby', np.float64), " <<
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
  w3 = 1.0;

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

void SlottedDiscs::printConstraints(void) const
{
  cout.setf(ios::scientific, ios::floatfield);
  cout.setf(ios::adjustfield, ios::right);
  cout.precision(9);
  for (int i = 0; i < 3; ++i) {
    cout.width(18);
    cout << con[i];
  }
  cout << endl;
} // printConstraints()

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
  z[121] = g*m;
  z[7] = cos(alpha);
  z[8] = sin(alpha);
} //evalConstants()

void SlottedDiscs::eoms(void)
{
  z[5] = cos(q3);
  z[6] = sin(q3);
  z[3] = cos(q2);
  z[11] = z[3]*z[6];
  z[4] = sin(q2);
  z[13] = 1 - pow(z[4],2);
  z[14] = pow(z[13],0.5);
  z[15] = 1/z[14];
  z[17] = ra*z[11]*z[15];
  z[12] = z[3]*z[5];
  z[19] = ra*z[12]*z[15];
  z[25] = z[7]*z[4] + z[8]*z[11];
  z[26] = 1 - pow(z[25],2);
  z[27] = pow(z[26],0.5);
  z[28] = 1/z[27];
  z[32] = rb*z[12]*z[28];
  z[24] = z[8]*z[4] - z[7]*z[11];
  z[30] = rb*z[24]*z[28];
  z[29] = z[25]/z[27];
  z[31] = rb*(z[29]-z[25]*z[28]);
  z[34] = -z[6]*z[17] - z[5]*(l+z[19]-z[32]) - z[6]*(z[7]*z[30]+z[8]*z[31]);
  z[16] = z[4]/z[14];
  z[18] = ra*(z[16]-z[4]*z[15]);
  z[38] = -z[4]*z[17] - z[11]*z[18] - z[24]*z[31] - z[25]*z[30];
  z[36] = z[11]*(l+z[19]-z[32]) - z[12]*z[17] - z[12]*(z[7]*z[30]+z[8]*z[31]);
  z[21] = z[8]*z[5];
  z[20] = z[7]*z[5];
  z[37] = z[5]*z[18] + z[21]*z[30] - z[20]*z[31];
  z[33] = z[6]*(z[7]*z[31]-z[18]-z[8]*z[30]);
  z[35] = z[12]*z[18] + z[4]*(l+z[19]-z[32]) - z[12]*(z[7]*z[31]-z[8]*z[30]);
  z[39] = -z[33]*z[36] - z[34]*z[35];
  z[40] = (z[34]*z[38]-z[36]*z[37])/z[39];
  z[42] = (z[5]+z[6]*z[40])/z[3];
  q1p = w3*z[42];
  z[43] = z[6] - z[5]*z[40];
  q2p = w3*z[43];
  z[10] = z[4]*z[5];
  z[9] = z[4]*z[6];
  z[23] = z[7]*z[3] - z[8]*z[9];
  z[53] = z[23]*z[25];
  z[55] = z[53]/pow(z[26],0.5);
  z[57] = z[55]/pow(z[27],2);
  z[65] = rb*(z[10]*z[28]-z[12]*z[57]);
  z[61] = (z[23]*z[27]+z[25]*z[55])/pow(z[27],2);
  z[63] = rb*(z[61]-z[23]*z[28]-z[25]*z[57]);
  z[44] = z[3]*z[4];
  z[45] = z[44]/pow(z[13],0.5);
  z[48] = (z[3]*z[14]+z[4]*z[45])/pow(z[14],2);
  z[46] = z[45]/pow(z[14],2);
  z[49] = ra*(z[48]-z[3]*z[15]-z[4]*z[46]);
  z[47] = ra*(z[9]*z[15]-z[11]*z[46]);
  z[50] = ra*(z[10]*z[15]-z[12]*z[46]);
  z[22] = z[7]*z[9] + z[8]*z[3];
  z[59] = rb*(z[22]*z[28]+z[24]*z[57]);
  z[68] = z[10]*z[32] + z[12]*z[65] + z[23]*z[31] + z[25]*z[63] - z[3]*z[18] - 
  z[4]*z[49] - z[9]*z[17] - z[11]*z[47] - z[12]*z[50] - z[22]*z[30] - z[24]*
  z[59] - z[10]*(l+z[19]);
  z[52] = z[8]*z[12];
  z[54] = z[25]*z[52];
  z[56] = z[54]/pow(z[26],0.5);
  z[58] = z[56]/pow(z[27],2);
  z[51] = z[7]*z[12];
  z[60] = rb*(z[24]*z[58]-z[28]*z[51]);
  z[66] = rb*(z[11]*z[28]-z[12]*z[58]);
  z[62] = (z[25]*z[56]+z[27]*z[52])/pow(z[27],2);
  z[64] = rb*(z[62]-z[25]*z[58]-z[28]*z[52]);
  z[67] = l*z[11] + z[24]*z[60] - z[11]*z[32] - z[12]*z[66] - z[25]*z[64] - 
  z[30]*z[51] - z[31]*z[52];
  z[69] = z[68]/z[67];
  z[70] = z[43]*z[69];
  q3p = w3*z[70];
  z[1] = cos(q1);
  z[2] = sin(q1);
  q4p = -(z[17]*(z[1]*z[6]+z[2]*z[10])+z[19]*(z[1]*z[5]-z[2]*z[9]))*q3p;
  q5p = -(z[19]*(z[1]*z[9]+z[2]*z[5])-z[17]*(z[1]*z[10]-z[2]*z[6]))*q3p;
  z[41] = (z[33]*z[38]+z[35]*z[37])/z[39];
  z[71] = z[41]*(k+z[19]) - z[18];
  z[72] = z[17] - z[40]*(k+z[19]);
  z[73] = z[17]*z[41] - z[18]*z[40];
  z[80] = z[10]*z[17] + z[12]*z[47] + z[10]*(z[7]*z[30]+z[8]*z[31]) - z[11]*(
  z[50]-z[65]) - z[9]*(l+z[19]-z[32]) - z[12]*(z[7]*z[59]+z[8]*z[63]);
  z[88] = z[12]*z[49] + z[3]*(l+z[19]-z[32]) + z[10]*(z[7]*z[31]-z[8]*z[30]) - 
  z[10]*z[18] - z[4]*(z[50]-z[65]) - z[12]*(z[7]*z[63]-z[8]*z[59]);
  z[77] = z[6]*z[47] + z[5]*(z[50]-z[65]) - z[6]*(z[7]*z[59]+z[8]*z[63]);
  z[87] = z[6]*(z[7]*z[63]-z[49]-z[8]*z[59]);
  z[90] = -z[33]*z[80] - z[34]*z[88] - z[35]*z[77] - z[36]*z[87];
  z[78] = z[4]*z[47] + z[9]*z[18] - z[3]*z[17] - z[11]*z[49] - z[22]*z[31] - 
  z[23]*z[30] - z[24]*z[63] - z[25]*z[59];
  z[84] = z[5]*z[49] + z[21]*z[59] - z[20]*z[63];
  z[96] = (z[90]*(z[33]*z[38]+z[35]*z[37])-z[39]*(z[33]*z[78]+z[35]*z[84]+
  z[37]*z[88]+z[38]*z[87]))/pow(z[39],2);
  z[105] = -z[49] - z[41]*z[50] - z[96]*(k+z[19]);
  z[107] = w3*z[105];
  z[81] = z[11]*z[66] + z[12]*(l-z[32]) + z[11]*(z[7]*z[30]+z[8]*z[31]) - 
  z[12]*(z[7]*z[60]+z[8]*z[64]);
  z[89] = z[11]*(z[7]*z[31]-z[8]*z[30]) - z[11]*z[18] - z[4]*(z[17]-z[66]) - 
  z[12]*(z[7]*z[64]-z[8]*z[60]);
  z[76] = z[6]*(l-z[32]) - z[5]*z[66] - z[5]*(z[7]*z[30]+z[8]*z[31]) - z[6]*(
  z[7]*z[60]+z[8]*z[64]);
  z[86] = z[6]*(z[7]*z[64]-z[8]*z[60]) + z[5]*(z[7]*z[31]-z[18]-z[8]*z[30]);
  z[91] = -z[33]*z[81] - z[34]*z[89] - z[35]*z[76] - z[36]*z[86];
  z[79] = z[31]*z[51] - z[4]*z[19] - z[12]*z[18] - z[24]*z[64] - z[25]*z[60] - 
  z[30]*z[52];
  z[83] = z[7]*z[6];
  z[82] = z[8]*z[6];
  z[85] = z[21]*z[60] + z[31]*z[83] - z[6]*z[18] - z[20]*z[64] - z[30]*z[82];
  z[97] = (z[91]*(z[33]*z[38]+z[35]*z[37])-z[39]*(z[33]*z[79]+z[35]*z[85]+
  z[37]*z[89]+z[38]*z[86]))/pow(z[39],2);
  z[106] = -z[17]*z[41] - z[97]*(k+z[19]);
  z[108] = w3*z[106];
  z[117] = z[107]*q2p + z[108]*q3p - pow(w3,2)*z[72] - pow(w3,2)*z[41]*z[73];
  z[93] = (z[91]*(z[34]*z[38]-z[36]*z[37])+z[39]*(z[36]*z[85]+z[37]*z[81]-
  z[34]*z[79]-z[38]*z[76]))/pow(z[39],2);
  z[109] = z[19] + z[17]*z[40] + z[93]*(k+z[19]);
  z[111] = w3*z[109];
  z[92] = (z[90]*(z[34]*z[38]-z[36]*z[37])+z[39]*(z[36]*z[84]+z[37]*z[80]-
  z[34]*z[78]-z[38]*z[77]))/pow(z[39],2);
  z[110] = z[40]*z[50] + z[92]*(k+z[19]) - z[47];
  z[112] = w3*z[110];
  z[118] = pow(w3,2)*z[71] + pow(w3,2)*z[40]*z[73] + z[111]*q3p + z[112]*q2p;
  z[113] = z[18]*z[92] - z[17]*z[96] - z[40]*z[49] - z[41]*z[47];
  z[115] = w3*z[113];
  z[114] = z[18]*z[93] + z[19]*z[41] - z[17]*z[97];
  z[116] = w3*z[114];
  z[119] = pow(w3,2)*z[41]*z[71] + z[115]*q2p + z[116]*q3p - pow(w3,2)*z[40]*
  z[72];
  z[94] = w3*z[92];
  z[95] = w3*z[93];
  z[100] = z[94]*q2p + z[95]*q3p;
  z[98] = w3*z[96];
  z[99] = w3*z[97];
  z[101] = z[98]*q2p + z[99]*q3p;
  z[122] = z[121]*(z[11]*z[71]-z[4]*z[72]-z[12]*z[73]) + m*(z[71]*z[117]+
  z[72]*z[118]+z[73]*z[119]) - z[40]*(Ixx*z[100]+Ixy*z[101]) - z[41]*(Ixy*
  z[100]+Iyy*z[101]);
  z[120] = -Izz - Iyy*pow(z[41],2) - z[40]*(Ixx*z[40]+2*Ixy*z[41]) - m*(pow(
  z[71],2)+pow(z[72],2)+pow(z[73],2));
  z[123] = z[122]/z[120];
  w3p = z[123];
} // eoms()

void SlottedDiscs::computeOutputs(void)
{
  w1 = -w3*z[40];
  w2 = -w3*z[41];
  ke = 0.5*pow(w3,2)*(Izz+Iyy*pow(z[41],2)+z[40]*(Ixx*z[40]+2*Ixy*z[41])+m*(
  pow(z[71],2)+pow(z[72],2)+pow(z[73],2)));
  pe = -g*m*(z[4]*z[18]-z[11]*z[17]-z[12]*(k+z[19]));
  te = ke + pe;
  z[124] = z[1]*z[5] - z[2]*z[9];
  z[125] = z[1]*z[9] + z[2]*z[5];
  z[126] = z[2]*z[3];
  z[127] = z[1]*z[3];
  z[128] = z[1]*z[6] + z[2]*z[10];
  z[129] = z[2]*z[6] - z[1]*z[10];
  z[130] = z[7]*z[124] - z[8]*z[126];
  z[131] = z[7]*z[125] + z[8]*z[127];
  z[132] = -z[7]*z[126] - z[8]*z[124];
  z[133] = z[7]*z[127] - z[8]*z[125];

  con[0] = w3*z[20]*z[31] + z[6]*(w1*z[18]-w2*z[17]) + z[5]*(w2*z[32]-w3*
  z[18]-w2*(l+z[19])) + z[6]*(z[8]*w1*z[30]-z[7]*w1*z[31]-z[7]*w2*z[30]-z[8]*
  w2*z[31]) - w3*z[21]*z[30];
  con[1] = w3*z[24]*z[31] + w3*z[25]*z[30] + z[12]*(w1*z[18]-w2*z[17]) + 
  z[12]*(z[8]*w1*z[30]-z[7]*w1*z[31]-z[7]*w2*z[30]-z[8]*w2*z[31]) - z[4]*(w1*
  z[32]-w3*z[17]-w1*(l+z[19])) - z[11]*(w2*z[32]-w3*z[18]-w2*(l+z[19]));
  con[2] = z[4]*z[18] + z[12]*z[32] + z[24]*z[30] - z[11]*z[17] - z[25]*z[31] - 
  z[12]*(l+z[19]);
  no_cb[0] = q4 + z[17]*z[124] + z[30]*z[130] + z[32]*z[128] - z[18]*z[126] - 
  z[31]*z[132] - z[128]*(l+z[19]);
  no_cb[1] = q5 + z[17]*z[125] + z[18]*z[127] + z[30]*z[131] + z[32]*z[129] - 
  z[31]*z[133] - z[129]*(l+z[19]);
} // computeOutputs()

void SlottedDiscs::setState(const double state[6])
{
  q1 = state[0];
  q2 = state[1];
  q3 = state[2];
  q4 = state[3];
  q5 = state[4];
  w3 = state[5];
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
