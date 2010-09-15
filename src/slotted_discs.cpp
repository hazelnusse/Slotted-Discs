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
  z[109] = g*m;
  z[7] = cos(alpha);
  z[8] = sin(alpha);
} //evalConstants()

void SlottedDiscs::eoms(void)
{
  z[5] = cos(q3);
  z[6] = sin(q3);
  z[9] = ra*z[6];
  z[10] = ra*z[5];
  z[3] = cos(q2);
  z[14] = z[3]*z[5];
  z[4] = sin(q2);
  z[13] = z[3]*z[6];
  z[20] = z[7]*z[4] + z[8]*z[13];
  z[21] = 1 - pow(z[20],2);
  z[22] = pow(z[21],0.5);
  z[23] = 1/z[22];
  z[27] = rb*z[14]*z[23];
  z[19] = z[8]*z[4] - z[7]*z[13];
  z[25] = rb*z[19]*z[23];
  z[24] = z[20]/z[22];
  z[26] = rb*(z[24]-z[20]*z[23]);
  z[29] = -z[6]*z[9] - z[5]*(l+z[10]-z[27]) - z[6]*(z[7]*z[25]+z[8]*z[26]);
  z[33] = -z[4]*z[9] - z[19]*z[26] - z[20]*z[25];
  z[31] = z[13]*(l+z[10]-z[27]) - z[9]*z[14] - z[14]*(z[7]*z[25]+z[8]*z[26]);
  z[16] = z[8]*z[5];
  z[15] = z[7]*z[5];
  z[32] = z[16]*z[25] - z[15]*z[26];
  z[28] = z[6]*(z[7]*z[26]-z[8]*z[25]);
  z[30] = z[4]*(l+z[10]-z[27]) - z[14]*(z[7]*z[26]-z[8]*z[25]);
  z[34] = -z[28]*z[31] - z[29]*z[30];
  z[35] = (z[29]*z[33]-z[31]*z[32])/z[34];
  z[37] = (z[5]+z[6]*z[35])/z[3];
  q1p = w3*z[37];
  z[38] = z[6] - z[5]*z[35];
  q2p = w3*z[38];
  z[12] = z[4]*z[5];
  z[11] = z[4]*z[6];
  z[18] = z[7]*z[3] - z[8]*z[11];
  z[41] = z[18]*z[20];
  z[43] = z[41]/pow(z[21],0.5);
  z[45] = z[43]/pow(z[22],2);
  z[53] = rb*(z[12]*z[23]-z[14]*z[45]);
  z[49] = (z[18]*z[22]+z[20]*z[43])/pow(z[22],2);
  z[51] = rb*(z[49]-z[18]*z[23]-z[20]*z[45]);
  z[17] = z[7]*z[11] + z[8]*z[3];
  z[47] = rb*(z[17]*z[23]+z[19]*z[45]);
  z[56] = z[12]*z[27] + z[14]*z[53] + z[18]*z[26] + z[20]*z[51] - z[9]*z[11] - 
  z[17]*z[25] - z[19]*z[47] - z[12]*(l+z[10]);
  z[40] = z[8]*z[14];
  z[42] = z[20]*z[40];
  z[44] = z[42]/pow(z[21],0.5);
  z[46] = z[44]/pow(z[22],2);
  z[39] = z[7]*z[14];
  z[48] = rb*(z[19]*z[46]-z[23]*z[39]);
  z[54] = rb*(z[13]*z[23]-z[14]*z[46]);
  z[50] = (z[20]*z[44]+z[22]*z[40])/pow(z[22],2);
  z[52] = rb*(z[50]-z[20]*z[46]-z[23]*z[40]);
  z[55] = l*z[13] + z[19]*z[48] - z[13]*z[27] - z[14]*z[54] - z[20]*z[52] - 
  z[25]*z[39] - z[26]*z[40];
  z[57] = z[56]/z[55];
  z[58] = z[38]*z[57];
  q3p = w3*z[58];
  z[1] = cos(q1);
  z[2] = sin(q1);
  q4p = -(z[9]*(z[1]*z[6]+z[2]*z[12])+z[10]*(z[1]*z[5]-z[2]*z[11]))*q3p;
  q5p = -(z[10]*(z[1]*z[11]+z[2]*z[5])-z[9]*(z[1]*z[12]-z[2]*z[6]))*q3p;
  z[36] = (z[28]*z[33]+z[30]*z[32])/z[34];
  z[59] = z[36]*(k+z[10]);
  z[60] = z[9] - z[35]*(k+z[10]);
  z[61] = z[9]*z[36];
  z[69] = z[13]*z[54] + z[14]*(l-z[27]) + z[13]*(z[7]*z[25]+z[8]*z[26]) - 
  z[14]*(z[7]*z[48]+z[8]*z[52]);
  z[77] = z[13]*(z[7]*z[26]-z[8]*z[25]) - z[4]*(z[9]-z[54]) - z[14]*(z[7]*
  z[52]-z[8]*z[48]);
  z[64] = z[6]*(l-z[27]) - z[5]*z[54] - z[5]*(z[7]*z[25]+z[8]*z[26]) - z[6]*(
  z[7]*z[48]+z[8]*z[52]);
  z[74] = z[5]*(z[7]*z[26]-z[8]*z[25]) + z[6]*(z[7]*z[52]-z[8]*z[48]);
  z[79] = -z[28]*z[69] - z[29]*z[77] - z[30]*z[64] - z[31]*z[74];
  z[66] = z[26]*z[39] - z[4]*z[10] - z[19]*z[52] - z[20]*z[48] - z[25]*z[40];
  z[71] = z[7]*z[6];
  z[70] = z[8]*z[6];
  z[72] = z[16]*z[48] + z[26]*z[71] - z[15]*z[52] - z[25]*z[70];
  z[85] = (z[79]*(z[28]*z[33]+z[30]*z[32])-z[34]*(z[28]*z[66]+z[30]*z[72]+
  z[32]*z[77]+z[33]*z[74]))/pow(z[34],2);
  z[93] = -z[9]*z[36] - z[85]*(k+z[10]);
  z[95] = w3*z[93];
  z[68] = z[9]*z[12] + z[13]*z[53] + z[12]*(z[7]*z[25]+z[8]*z[26]) - z[11]*(l+
  z[10]-z[27]) - z[14]*(z[7]*z[47]+z[8]*z[51]);
  z[76] = z[4]*z[53] + z[3]*(l+z[10]-z[27]) + z[12]*(z[7]*z[26]-z[8]*z[25]) - 
  z[14]*(z[7]*z[51]-z[8]*z[47]);
  z[65] = -z[5]*z[53] - z[6]*(z[7]*z[47]+z[8]*z[51]);
  z[75] = z[6]*(z[7]*z[51]-z[8]*z[47]);
  z[78] = -z[28]*z[68] - z[29]*z[76] - z[30]*z[65] - z[31]*z[75];
  z[67] = -z[3]*z[9] - z[17]*z[26] - z[18]*z[25] - z[19]*z[51] - z[20]*z[47];
  z[73] = z[16]*z[47] - z[15]*z[51];
  z[84] = (z[78]*(z[28]*z[33]+z[30]*z[32])-z[34]*(z[28]*z[67]+z[30]*z[73]+
  z[32]*z[76]+z[33]*z[75]))/pow(z[34],2);
  z[94] = z[84]*(k+z[10]);
  z[96] = w3*z[94];
  z[105] = z[95]*q3p - pow(w3,2)*z[60] - pow(w3,2)*z[36]*z[61] - z[96]*q2p;
  z[81] = (z[79]*(z[29]*z[33]-z[31]*z[32])+z[34]*(z[31]*z[72]+z[32]*z[69]-
  z[29]*z[66]-z[33]*z[64]))/pow(z[34],2);
  z[97] = z[10] + z[9]*z[35] + z[81]*(k+z[10]);
  z[99] = w3*z[97];
  z[80] = (z[78]*(z[29]*z[33]-z[31]*z[32])+z[34]*(z[31]*z[73]+z[32]*z[68]-
  z[29]*z[67]-z[33]*z[65]))/pow(z[34],2);
  z[98] = z[80]*(k+z[10]);
  z[100] = w3*z[98];
  z[106] = pow(w3,2)*z[59] + pow(w3,2)*z[35]*z[61] + z[99]*q3p + z[100]*q2p;
  z[101] = z[10]*z[36] - z[9]*z[85];
  z[103] = w3*z[101];
  z[102] = z[9]*z[84];
  z[104] = w3*z[102];
  z[107] = pow(w3,2)*z[36]*z[59] + z[103]*q3p - pow(w3,2)*z[35]*z[60] - 
  z[104]*q2p;
  z[82] = w3*z[80];
  z[83] = w3*z[81];
  z[88] = z[82]*q2p + z[83]*q3p;
  z[86] = w3*z[84];
  z[87] = w3*z[85];
  z[89] = z[86]*q2p + z[87]*q3p;
  z[110] = z[109]*(z[13]*z[59]-z[4]*z[60]-z[14]*z[61]) + m*(z[59]*z[105]+
  z[60]*z[106]+z[61]*z[107]) - z[35]*(Ixx*z[88]+Ixy*z[89]) - z[36]*(Ixy*z[88]+
  Iyy*z[89]);
  z[108] = -Izz - Iyy*pow(z[36],2) - z[35]*(Ixx*z[35]+2*Ixy*z[36]) - m*(pow(
  z[59],2)+pow(z[60],2)+pow(z[61],2));
  z[111] = z[110]/z[108];
  w3p = z[111];
} // eoms()

void SlottedDiscs::computeOutputs(void)
{
  w1 = -w3*z[35];
  w2 = -w3*z[36];
  ke = 0.5*pow(w3,2)*(Izz+Iyy*pow(z[36],2)+z[35]*(Ixx*z[35]+2*Ixy*z[36])+m*(
  pow(z[59],2)+pow(z[60],2)+pow(z[61],2)));
  pe = g*m*(z[9]*z[13]+z[14]*(k+z[10]));
  te = ke + pe;
  z[112] = z[1]*z[5] - z[2]*z[11];
  z[113] = z[1]*z[11] + z[2]*z[5];
  z[114] = z[2]*z[3];
  z[115] = z[1]*z[3];
  z[116] = z[1]*z[6] + z[2]*z[12];
  z[117] = z[2]*z[6] - z[1]*z[12];
  z[118] = z[7]*z[112] - z[8]*z[114];
  z[119] = z[7]*z[113] + z[8]*z[115];
  z[120] = -z[7]*z[114] - z[8]*z[112];
  z[121] = z[7]*z[115] - z[8]*z[113];

  con[0] = w3*z[15]*z[26] + z[6]*(z[8]*w1*z[25]-z[7]*w1*z[26]-z[7]*w2*z[25]-
  z[8]*w2*z[26]) - w2*z[6]*z[9] - w3*z[16]*z[25] - w2*z[5]*(l+z[10]-z[27]);
  con[1] = w3*z[19]*z[26] + w3*z[20]*z[25] + w2*z[13]*(l+z[10]-z[27]) + z[14]*(
  z[8]*w1*z[25]-z[7]*w1*z[26]-z[7]*w2*z[25]-z[8]*w2*z[26]) - w2*z[9]*z[14] - 
  z[4]*(w1*z[27]-w3*z[9]-w1*(l+z[10]));
  con[2] = z[14]*z[27] + z[19]*z[25] - z[9]*z[13] - z[20]*z[26] - z[14]*(l+
  z[10]);
  no_cb[0] = q4 + z[9]*z[112] + z[25]*z[118] + z[27]*z[116] - z[26]*z[120] - 
  z[116]*(l+z[10]);
  no_cb[1] = q5 + z[9]*z[113] + z[25]*z[119] + z[27]*z[117] - z[26]*z[121] - 
  z[117]*(l+z[10]);
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
