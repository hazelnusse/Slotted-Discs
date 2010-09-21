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
  file.write((char *) &(discs->H), sizeof discs->H);
  file.write((char *) &(discs->p), sizeof discs->p);
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
        "('te', np.float64), " <<
        "('H1', np.float64), " <<
        "('H2', np.float64), " <<
        "('H3', np.float64), " <<
        "('p1', np.float64), " <<
        "('p2', np.float64), " <<
        "('p3', np.float64)]) " << endl;
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
  z[108] = g*m;
  z[7] = cos(alpha);
  z[8] = sin(alpha);
} //evalConstants()

void SlottedDiscs::eoms(void)
{
  z[5] = cos(q3);
  z[3] = cos(q2);
  z[14] = z[3]*z[5];
  z[4] = sin(q2);
  z[6] = sin(q3);
  z[13] = z[3]*z[6];
  z[20] = z[7]*z[4] + z[8]*z[13];
  z[21] = 1 - pow(z[20],2);
  z[22] = pow(z[21],0.5);
  z[23] = 1/z[22];
  z[27] = rb*z[14]*z[23];
  z[10] = ra*z[5];
  z[29] = z[27] - l - z[10];
  z[24] = z[20]/z[22];
  z[26] = rb*(z[7]*z[24]-z[4]*z[23]);
  z[9] = ra*z[6];
  z[25] = rb*(z[8]*z[24]-z[13]*z[23]);
  z[28] = z[9] + z[25];
  z[30] = pow(z[26],2) + pow(z[28],2) + pow(z[29],2);
  z[31] = pow(z[30],0.5);
  z[33] = z[29]/z[31];
  z[32] = z[28]/z[31];
  z[35] = (z[5]*z[33]-z[6]*z[32])/z[3];
  q1p = w*z[35];
  z[36] = z[5]*z[32] + z[6]*z[33];
  q2p = w*z[36];
  z[34] = z[26]/z[31];
  z[37] = tan(q2);
  z[38] = -z[34] - z[37]*(z[5]*z[33]-z[6]*z[32]);
  q3p = w*z[38];
  z[1] = cos(q1);
  z[2] = sin(q1);
  z[12] = z[4]*z[5];
  z[43] = z[1]*z[6] + z[2]*z[12];
  z[11] = z[4]*z[6];
  z[39] = z[1]*z[5] - z[2]*z[11];
  z[45] = (z[9]*z[43]+z[10]*z[39])*q3p;
  q4p = -z[45];
  z[44] = z[2]*z[6] - z[1]*z[12];
  z[40] = z[1]*z[11] + z[2]*z[5];
  z[46] = (z[9]*z[44]+z[10]*z[40])*q3p;
  q5p = -z[46];
  z[47] = z[34]*(k+z[10]);
  z[48] = z[9]*z[33] + z[32]*(k+z[10]);
  z[49] = z[9]*z[34];
  z[18] = z[7]*z[3] - z[8]*z[11];
  z[54] = z[18]*z[20];
  z[56] = z[54]/pow(z[21],0.5);
  z[58] = (z[18]*z[22]+z[20]*z[56])/pow(z[22],2);
  z[60] = z[56]/pow(z[22],2);
  z[65] = rb*(z[7]*z[58]-z[3]*z[23]-z[4]*z[60]);
  z[62] = rb*(z[8]*z[58]+z[11]*z[23]-z[13]*z[60]);
  z[67] = rb*(z[12]*z[23]-z[14]*z[60]);
  z[70] = 2*z[26]*z[65] + 2*z[28]*z[62] - 2*z[29]*z[67];
  z[72] = z[70]/pow(z[30],0.5);
  z[82] = (z[29]*z[72]+2*z[31]*z[67])/pow(z[31],2);
  z[84] = w*z[82];
  z[53] = z[8]*z[14];
  z[55] = z[20]*z[53];
  z[57] = z[55]/pow(z[21],0.5);
  z[59] = (z[20]*z[57]+z[22]*z[53])/pow(z[22],2);
  z[61] = z[57]/pow(z[22],2);
  z[66] = rb*(z[7]*z[59]-z[4]*z[61]);
  z[63] = rb*(z[8]*z[59]-z[13]*z[61]-z[14]*z[23]);
  z[64] = z[10] + z[63];
  z[68] = rb*(z[13]*z[23]-z[14]*z[61]);
  z[69] = z[9] - z[68];
  z[71] = 2*z[26]*z[66] + 2*z[28]*z[64] + 2*z[29]*z[69];
  z[73] = z[71]/pow(z[30],0.5);
  z[83] = (z[29]*z[73]-2*z[31]*z[69])/pow(z[31],2);
  z[85] = w*z[83];
  z[88] = -0.5*z[84]*q2p - 0.5*z[85]*q3p;
  z[74] = (z[28]*z[72]-2*z[31]*z[62])/pow(z[31],2);
  z[76] = w*z[74];
  z[75] = (z[28]*z[73]-2*z[31]*z[64])/pow(z[31],2);
  z[77] = w*z[75];
  z[86] = -0.5*z[76]*q2p - 0.5*z[77]*q3p;
  z[78] = (z[26]*z[72]-2*z[31]*z[65])/pow(z[31],2);
  z[80] = w*z[78];
  z[79] = (z[26]*z[73]-2*z[31]*z[66])/pow(z[31],2);
  z[81] = w*z[79];
  z[87] = 0.5*z[80]*q2p + 0.5*z[81]*q3p;
  z[92] = -z[9]*z[34] - 0.5*z[79]*(k+z[10]);
  z[94] = w*z[92];
  z[93] = z[78]*(k+z[10]);
  z[95] = w*z[93];
  z[104] = z[94]*q3p - pow(w,2)*z[33]*z[48] - pow(w,2)*z[34]*z[49] - 0.5*
  z[95]*q2p;
  z[96] = z[10]*z[33] - z[9]*z[32] - 0.5*z[9]*z[83] - 0.5*z[75]*(k+z[10]);
  z[98] = w*z[96];
  z[97] = -0.5*z[9]*z[82] - 0.5*z[74]*(k+z[10]);
  z[99] = w*z[97];
  z[105] = pow(w,2)*z[33]*z[47] + z[98]*q3p + z[99]*q2p - pow(w,2)*z[32]*
  z[49];
  z[100] = z[10]*z[34] - 0.5*z[9]*z[79];
  z[102] = w*z[100];
  z[101] = z[9]*z[78];
  z[103] = w*z[101];
  z[106] = pow(w,2)*z[32]*z[48] + pow(w,2)*z[34]*z[47] + z[102]*q3p - 0.5*
  z[103]*q2p;
  z[109] = z[108]*(z[13]*z[47]-z[4]*z[48]-z[14]*z[49]) + Izz*z[33]*z[88] + 
  z[32]*(Ixx*z[86]+Ixy*z[87]) + m*(z[47]*z[104]+z[48]*z[105]+z[49]*z[106]) - 
  z[34]*(Ixy*z[86]+Iyy*z[87]);
  z[107] = -Iyy*pow(z[34],2) - Izz*pow(z[33],2) - z[32]*(Ixx*z[32]-2*Ixy*
  z[34]) - m*(pow(z[47],2)+pow(z[48],2)+pow(z[49],2));
  z[110] = z[109]/z[107];
  wp = z[110];
} // eoms()

void SlottedDiscs::computeOutputs(void)
{
  z[41] = z[2]*z[3];
  z[42] = z[1]*z[3];
  ke = 0.5*pow(w,2)*(Iyy*pow(z[34],2)+Izz*pow(z[33],2)+z[32]*(Ixx*z[32]-2*Ixy*
  z[34])+m*(pow(z[47],2)+pow(z[48],2)+pow(z[49],2)));
  pe = g*m*(z[9]*z[13]+z[14]*(k+z[10]));
  te = ke + pe;
  w1 = w*z[32];
  w2 = -w*z[34];
  w3 = w*z[33];
  z[111] = z[7]*z[3] - z[8]*z[4]*z[6];
  z[112] = z[20]*z[111]/(pow(z[21],0.5)*pow(z[22],2));
  z[113] = rb*(z[14]*z[112]-z[4]*z[5]*z[23]);
  z[114] = z[111]*(z[22]+pow(z[20],2)/pow(z[21],0.5))/pow(z[22],2);
  z[115] = rb*(z[7]*z[114]-z[3]*z[23]-z[4]*z[112]);
  z[116] = rb*(z[13]*z[112]-z[8]*z[114]-z[4]*z[6]*z[23]);
  z[117] = 2*z[26]*z[115] + 2*z[29]*z[113] - 2*z[28]*z[116];
  z[118] = (2*z[31]*z[113]-z[29]*z[117]/pow(z[30],0.5))/pow(z[31],2);
  z[119] = (2*z[31]*z[116]+z[28]*z[117]/pow(z[30],0.5))/pow(z[31],2);
  z[120] = (z[3]*(z[5]*z[118]+z[6]*z[119])+2*z[4]*(z[5]*z[33]-z[6]*z[32]))/
  pow(z[3],2);
  z[121] = z[8]*z[3]*z[5]*z[20];
  z[122] = rb*(z[3]*z[6]*z[23]-z[14]*z[121]/(pow(z[21],0.5)*pow(z[22],2)));
  z[123] = (z[8]*z[3]*z[5]*z[22]+z[20]*z[121]/pow(z[21],0.5))/pow(z[22],2);
  z[124] = rb*(z[7]*z[123]-z[4]*z[121]/(pow(z[21],0.5)*pow(z[22],2)));
  z[125] = rb*(z[8]*z[123]-z[3]*z[5]*z[23]-z[13]*z[121]/(pow(z[21],0.5)*pow(
  z[22],2)));
  z[126] = 2*z[26]*z[124] + 2*z[28]*(z[125]+ra*z[5]) - 2*z[29]*(z[122]-ra*
  z[6]);
  z[127] = (z[29]*z[126]/pow(z[30],0.5)+2*z[31]*(z[122]-ra*z[6]))/pow(z[31],2);
  z[128] = (z[28]*z[126]/pow(z[30],0.5)-2*z[31]*(z[125]+ra*z[5]))/pow(z[31],2);
  z[129] = (z[5]*z[127]+2*z[5]*z[32]+2*z[6]*z[33]-z[6]*z[128])/z[3];
  z[130] = 0.5*z[6]*z[118] - 0.5*z[5]*z[119];
  z[131] = z[5]*z[33] - z[6]*z[32] - 0.5*z[5]*z[128] - 0.5*z[6]*z[127];
  z[132] = (2*z[31]*z[115]-z[26]*z[117]/pow(z[30],0.5))/pow(z[31],2);
  z[133] = -0.5*z[132] - 0.5*z[37]*(z[5]*z[118]+z[6]*z[119]) - (z[5]*z[33]-
  z[6]*z[32])/pow(z[3],2);
  z[134] = (2*z[31]*z[124]-z[26]*z[126]/pow(z[30],0.5))/pow(z[31],2);
  z[135] = 0.5*z[37]*(z[5]*z[127]+2*z[5]*z[32]+2*z[6]*z[33]-z[6]*z[128]) - 
  0.5*z[134];
  z[136] = z[1]*z[12] - z[2]*z[6];
  z[137] = -z[1]*z[11] - z[2]*z[5];
  z[138] = (z[9]*z[136]+z[10]*z[137])*q3p;
  z[139] = w*z[133]*(z[9]*z[43]+z[10]*z[39]) + z[2]*z[3]*(z[5]*z[9]-z[6]*
  z[10])*q3p;
  z[140] = z[1]*z[5] - z[2]*z[4]*z[6];
  z[141] = -z[1]*z[6] - z[2]*z[4]*z[5];
  z[142] = w*z[135]*(z[9]*z[43]+z[10]*z[39]) + (z[9]*z[140]+z[10]*z[141]+ra*
  z[5]*z[43]-ra*z[6]*z[39])*q3p;
  z[143] = z[38]*(z[9]*z[43]+z[10]*z[39]);
  z[144] = w*z[133]*(z[9]*z[44]+z[10]*z[40]) - z[1]*z[3]*(z[5]*z[9]-z[6]*
  z[10])*q3p;
  z[145] = z[2]*z[5] + z[1]*z[4]*z[6];
  z[146] = z[1]*z[4]*z[5] - z[2]*z[6];
  z[147] = w*z[135]*(z[9]*z[44]+z[10]*z[40]) + (z[9]*z[145]+z[10]*z[146]+ra*
  z[5]*z[44]-ra*z[6]*z[40])*q3p;
  z[148] = z[38]*(z[9]*z[44]+z[10]*z[40]);
  z[149] = 0.5*z[9]*z[118] - 0.5*z[119]*(k+z[10]);
  z[150] = -z[7]*z[4] - z[8]*z[3]*z[6];
  z[151] = z[18]*z[111] + z[20]*z[150];
  z[152] = (z[21]*z[151]+z[20]*z[54]*z[111])/pow(z[21],1.5);
  z[153] = (z[22]*(z[20]*z[152]+z[22]*z[150]+z[56]*z[111])+z[20]*z[111]*(
  z[18]*z[22]+2*z[20]*z[56])/pow(z[21],0.5))/pow(z[22],3);
  z[154] = (z[22]*z[152]+2*z[20]*z[56]*z[111]/pow(z[21],0.5))/pow(z[22],3);
  z[155] = rb*(z[7]*z[153]+z[4]*z[23]-z[3]*z[60]-z[3]*z[112]-z[4]*z[154]);
  z[156] = rb*(z[13]*z[154]-z[8]*z[153]-z[11]*z[112]-z[3]*z[6]*z[23]-z[4]*
  z[6]*z[60]);
  z[157] = rb*(z[14]*z[154]-z[12]*z[112]-z[3]*z[5]*z[23]-z[4]*z[5]*z[60]);
  z[158] = 2*z[26]*z[155] + 2*z[29]*z[157] + 2*z[65]*z[115] - 2*z[28]*z[156] - 
  2*z[62]*z[116] - 2*z[67]*z[113];
  z[159] = (2*z[30]*z[158]-z[70]*z[117])/pow(z[30],1.5);
  z[160] = (2*z[117]*(z[29]*z[72]+z[31]*z[67])/pow(z[30],0.5)+z[31]*(4*z[31]*
  z[157]-2*z[72]*z[113]-z[29]*z[159]))/pow(z[31],3);
  z[161] = z[53]*z[111] - z[8]*z[4]*z[5]*z[20];
  z[162] = (z[21]*z[161]+z[20]*z[55]*z[111])/pow(z[21],1.5);
  z[163] = (z[20]*z[111]*(z[22]*z[53]+2*z[20]*z[57])/pow(z[21],0.5)+z[22]*(
  z[20]*z[162]+z[57]*z[111]-z[8]*z[4]*z[5]*z[22]))/pow(z[22],3);
  z[164] = (z[22]*z[162]+2*z[20]*z[57]*z[111]/pow(z[21],0.5))/pow(z[22],3);
  z[165] = rb*(z[7]*z[163]-z[3]*z[61]-z[4]*z[164]);
  z[166] = rb*(z[13]*z[164]+z[14]*z[112]-z[8]*z[163]-z[4]*z[5]*z[23]-z[4]*
  z[6]*z[61]);
  z[167] = rb*(z[13]*z[112]+z[4]*z[5]*z[61]-z[14]*z[164]-z[4]*z[6]*z[23]);
  z[168] = 2*z[26]*z[165] + 2*z[66]*z[115] + 2*z[69]*z[113] - 2*z[28]*z[166] - 
  2*z[29]*z[167] - 2*z[64]*z[116];
  z[169] = (2*z[30]*z[168]-z[71]*z[117])/pow(z[30],1.5);
  z[170] = (2*z[29]*z[73]*z[117]/pow(z[30],0.5)-z[31]*(z[29]*z[169]+2*z[73]*
  z[113]+4*z[31]*z[167]+2*z[69]*z[117]/pow(z[30],0.5)))/pow(z[31],3);
  z[171] = w*(2*z[84]*z[130]+2*z[85]*z[133]-z[160]*q2p-z[170]*q3p);
  z[172] = (2*z[28]*z[72]*z[117]/pow(z[30],0.5)+z[31]*(2*z[72]*z[116]-4*z[31]*
  z[156]-z[28]*z[159]-2*z[62]*z[117]/pow(z[30],0.5)))/pow(z[31],3);
  z[173] = (2*z[28]*z[73]*z[117]/pow(z[30],0.5)+z[31]*(2*z[73]*z[116]-4*z[31]*
  z[166]-z[28]*z[169]-2*z[64]*z[117]/pow(z[30],0.5)))/pow(z[31],3);
  z[174] = w*(2*z[76]*z[130]+2*z[77]*z[133]-z[172]*q2p-z[173]*q3p);
  z[175] = (2*z[26]*z[72]*z[117]/pow(z[30],0.5)+z[31]*(4*z[31]*z[155]-2*z[72]*
  z[115]-z[26]*z[159]-2*z[65]*z[117]/pow(z[30],0.5)))/pow(z[31],3);
  z[176] = (2*z[26]*z[73]*z[117]/pow(z[30],0.5)+z[31]*(4*z[31]*z[165]-2*z[73]*
  z[115]-z[26]*z[169]-2*z[66]*z[117]/pow(z[30],0.5)))/pow(z[31],3);
  z[177] = w*(2*z[80]*z[130]+2*z[81]*z[133]-z[175]*q2p-z[176]*q3p);
  z[178] = 0.25*z[176]*(k+z[10]) - 0.5*z[9]*z[132];
  z[179] = w*z[175]*(k+z[10]);
  z[180] = w*z[94]*z[133] + 0.25*z[179]*q2p + w*z[178]*q3p - 0.5*w*z[95]*
  z[130] - pow(w,2)*z[33]*z[149] - 0.5*pow(w,2)*z[48]*z[118] - 0.5*pow(w,2)*
  z[49]*z[132] - 0.5*pow(w,2)*z[9]*z[34]*z[132];
  z[181] = 0.25*z[9]*z[170] + 0.5*z[9]*z[119] + 0.5*z[10]*z[118] + 0.25*
  z[173]*(k+z[10]);
  z[182] = 0.25*z[9]*z[160] + 0.25*z[172]*(k+z[10]);
  z[183] = w*(w*z[9]*z[32]*z[132]-2*z[98]*z[133]-2*z[99]*z[130]-w*z[47]*
  z[118]-w*z[49]*z[119]-w*z[33]*z[132]*(k+z[10])-2*z[181]*q3p-2*z[182]*q2p);
  z[184] = 0.25*z[9]*z[176] + 0.5*z[10]*z[132];
  z[185] = w*(2*z[103]*z[130]+2*w*z[48]*z[119]-4*z[102]*z[133]-4*w*z[32]*
  z[149]-2*w*z[47]*z[132]-2*w*z[34]*z[132]*(k+z[10])-4*z[184]*q3p-z[9]*z[175]*
  q2p);
  z[186] = 0.5*Izz*z[118]*z[88] + 0.25*z[34]*(Ixy*z[174]-Iyy*z[177]) - 0.5*
  z[108]*(2*z[3]*z[48]+2*z[4]*z[149]+z[9]*z[14]*z[132]+2*z[4]*z[6]*z[47]-2*
  z[4]*z[5]*z[49]-z[13]*z[132]*(k+z[10])) - 0.25*Izz*z[33]*z[171] - 0.5*
  z[119]*(Ixx*z[86]+Ixy*z[87]) - 0.5*z[132]*(Ixy*z[86]+Iyy*z[87]) - 0.25*
  z[32]*(Ixx*z[174]-Ixy*z[177]) - 0.25*m*(z[49]*z[185]+2*z[48]*z[183]-4*z[47]*
  z[180]-4*z[149]*z[105]-2*z[9]*z[132]*z[106]-2*z[132]*(k+z[10])*z[104]);
  z[187] = z[32]*(Ixx*z[119]+Ixy*z[132]) - Ixy*z[34]*z[119] - Iyy*z[34]*
  z[132] - Izz*z[33]*z[118] - m*(2*z[48]*z[149]+z[9]*z[49]*z[132]+z[47]*
  z[132]*(k+z[10]));
  z[188] = (z[107]*z[186]-z[187]*z[109])/pow(z[107],2);
  z[189] = 0.5*z[134]*(k+z[10]) - ra*z[6]*z[34];
  z[190] = ra*z[5]*z[33] - 0.5*z[9]*z[127] - ra*z[6]*z[32] - 0.5*z[128]*(k+
  z[10]);
  z[191] = 0.5*z[9]*z[134] + ra*z[5]*z[34];
  z[192] = z[8]*z[5]*(z[3]*z[18]-z[4]*z[20]);
  z[193] = (z[21]*z[192]+z[54]*z[121])/pow(z[21],1.5);
  z[194] = (z[121]*(z[18]*z[22]+2*z[20]*z[56])/pow(z[21],0.5)+z[22]*(z[20]*
  z[193]+z[8]*z[3]*z[5]*z[56]-z[8]*z[4]*z[5]*z[22]))/pow(z[22],3);
  z[195] = (z[22]*z[193]+2*z[56]*z[121]/pow(z[21],0.5))/pow(z[22],3);
  z[196] = rb*(z[7]*z[194]-z[4]*z[195]-z[3]*z[121]/(pow(z[21],0.5)*pow(z[22],
  2)));
  z[197] = rb*(z[13]*z[195]+z[3]*z[5]*z[60]-z[8]*z[194]-z[4]*z[5]*z[23]-z[11]*
  z[121]/(pow(z[21],0.5)*pow(z[22],2)));
  z[198] = rb*(z[14]*z[195]+z[4]*z[6]*z[23]-z[3]*z[6]*z[60]-z[12]*z[121]/(
  pow(z[21],0.5)*pow(z[22],2)));
  z[199] = 2*z[26]*z[196] + 2*z[29]*z[198] + 2*z[65]*z[124] + 2*z[62]*(z[125]+
  ra*z[5]) + 2*z[67]*(z[122]-ra*z[6]) - 2*z[28]*z[197];
  z[200] = (2*z[30]*z[199]-z[70]*z[126])/pow(z[30],1.5);
  z[201] = (2*z[126]*(z[29]*z[72]+z[31]*z[67])/pow(z[30],0.5)-z[31]*(z[29]*
  z[200]-4*z[31]*z[198]-2*z[72]*(z[122]-ra*z[6])))/pow(z[31],3);
  z[202] = z[8]*z[3]*(z[5]*z[53]-z[6]*z[20]);
  z[203] = (z[21]*z[202]+z[55]*z[121])/pow(z[21],1.5);
  z[204] = (z[121]*(z[22]*z[53]+2*z[20]*z[57])/pow(z[21],0.5)+z[22]*(z[20]*
  z[203]+z[8]*z[3]*z[5]*z[57]-z[8]*z[3]*z[6]*z[22]))/pow(z[22],3);
  z[205] = (z[22]*z[203]+2*z[57]*z[121]/pow(z[21],0.5))/pow(z[22],3);
  z[206] = rb*(z[7]*z[204]-z[4]*z[205]);
  z[207] = rb*(z[8]*z[204]+z[3]*z[6]*z[23]-z[13]*z[205]-z[3]*z[5]*z[61]-z[14]*
  z[121]/(pow(z[21],0.5)*pow(z[22],2)));
  z[208] = rb*(z[14]*z[205]-z[3]*z[5]*z[23]-z[3]*z[6]*z[61]-z[13]*z[121]/(
  pow(z[21],0.5)*pow(z[22],2)));
  z[209] = 2*z[26]*z[206] + 2*z[66]*z[124] + 2*z[29]*(z[208]+ra*z[5]) + 2*
  z[64]*(z[125]+ra*z[5]) + 2*z[28]*(z[207]-ra*z[6]) - 2*z[69]*(z[122]-ra*z[6]);
  z[210] = (2*z[30]*z[209]-z[71]*z[126])/pow(z[30],1.5);
  z[211] = (2*z[126]*(z[29]*z[73]-z[31]*z[69])/pow(z[30],0.5)-z[31]*(z[29]*
  z[210]-4*z[31]*(z[208]+ra*z[5])-2*z[73]*(z[122]-ra*z[6])))/pow(z[31],3);
  z[212] = w*(2*z[84]*z[131]+2*z[85]*z[135]-z[201]*q2p-z[211]*q3p);
  z[213] = (2*z[126]*(z[28]*z[72]-z[31]*z[62])/pow(z[30],0.5)-z[31]*(z[28]*
  z[200]+4*z[31]*z[197]+2*z[72]*(z[125]+ra*z[5])))/pow(z[31],3);
  z[214] = (2*z[126]*(z[28]*z[73]-z[31]*z[64])/pow(z[30],0.5)-z[31]*(z[28]*
  z[210]+2*z[73]*(z[125]+ra*z[5])-4*z[31]*(z[207]-ra*z[6])))/pow(z[31],3);
  z[215] = w*(2*z[76]*z[131]+2*z[77]*z[135]-z[213]*q2p-z[214]*q3p);
  z[216] = (2*z[26]*z[72]*z[126]/pow(z[30],0.5)+z[31]*(4*z[31]*z[196]-2*z[72]*
  z[124]-z[26]*z[200]-2*z[65]*z[126]/pow(z[30],0.5)))/pow(z[31],3);
  z[217] = (2*z[26]*z[73]*z[126]/pow(z[30],0.5)+z[31]*(4*z[31]*z[206]-2*z[73]*
  z[124]-z[26]*z[210]-2*z[66]*z[126]/pow(z[30],0.5)))/pow(z[31],3);
  z[218] = w*(2*z[80]*z[131]+2*z[81]*z[135]-z[216]*q2p-z[217]*q3p);
  z[219] = 0.5*ra*z[6]*z[79] + 0.25*z[217]*(k+z[10]) - 0.5*z[9]*z[134] - ra*
  z[5]*z[34];
  z[220] = -ra*z[6]*z[78] - 0.5*z[216]*(k+z[10]);
  z[221] = w*(2*z[94]*z[135]+w*z[48]*z[127]+2*z[219]*q3p-z[95]*z[131]-2*w*
  z[33]*z[190]-2*w*z[34]*z[191]-w*z[49]*z[134]-z[220]*q2p);
  z[222] = 0.25*z[9]*z[211] + 0.5*z[9]*z[128] + 0.5*ra*z[6]*z[75] + 0.25*
  z[214]*(k+z[10]) - 0.5*z[10]*z[127] - ra*z[5]*z[32] - ra*z[6]*z[33] - 0.5*
  ra*z[5]*z[83];
  z[223] = 0.25*z[9]*z[201] + 0.5*ra*z[6]*z[74] + 0.25*z[213]*(k+z[10]) - 0.5*
  ra*z[5]*z[82];
  z[224] = w*(w*z[47]*z[127]+2*w*z[32]*z[191]-2*z[98]*z[135]-2*z[99]*z[131]-2*
  w*z[33]*z[189]-w*z[49]*z[128]-2*z[222]*q3p-2*z[223]*q2p);
  z[225] = 0.25*z[9]*z[217] + 0.5*z[10]*z[134] - ra*z[6]*z[34] - 0.5*ra*z[5]*
  z[79];
  z[226] = ra*z[5]*z[78] - 0.5*z[9]*z[216];
  z[227] = w*(2*z[102]*z[135]+w*z[47]*z[134]+2*w*z[32]*z[190]+2*w*z[34]*
  z[189]+2*z[225]*q3p-z[103]*z[131]-w*z[48]*z[128]-z[226]*q2p);
  z[228] = 0.25*z[34]*(Ixy*z[215]-Iyy*z[218]) - z[108]*(z[4]*z[190]+z[14]*
  z[191]-z[13]*z[189]-z[3]*z[5]*z[47]-z[3]*z[6]*z[49]) - 0.5*Izz*z[127]*z[88] - 
  0.25*Izz*z[33]*z[212] - 0.5*z[128]*(Ixx*z[86]+Ixy*z[87]) - 0.5*z[134]*(Ixy*
  z[86]+Iyy*z[87]) - 0.25*z[32]*(Ixx*z[215]-Ixy*z[218]) - 0.5*m*(z[48]*z[224]-
  2*z[189]*z[104]-2*z[190]*z[105]-2*z[191]*z[106]-z[47]*z[221]-z[49]*z[227]);
  z[229] = Izz*z[33]*z[127] + z[32]*(Ixx*z[128]+Ixy*z[134]) - Ixy*z[34]*
  z[128] - Iyy*z[34]*z[134] - 2*m*(z[47]*z[189]+z[48]*z[190]+z[49]*z[191]);
  z[230] = (z[107]*z[228]-z[229]*z[109])/pow(z[107],2);
  z[231] = -0.5*z[36]*z[84] - 0.5*z[38]*z[85] - 0.5*z[82]*q2p - 0.5*z[83]*q3p;
  z[232] = -0.5*z[36]*z[76] - 0.5*z[38]*z[77] - 0.5*z[74]*q2p - 0.5*z[75]*q3p;
  z[233] = 0.5*z[36]*z[80] + 0.5*z[38]*z[81] + 0.5*z[78]*q2p + 0.5*z[79]*q3p;
  z[234] = z[38]*z[94] + z[92]*q3p - 0.5*z[36]*z[95] - 2*w*z[33]*z[48] - 2*w*
  z[34]*z[49] - 0.5*z[93]*q2p;
  z[235] = z[36]*z[99] + z[38]*z[98] + 2*w*z[33]*z[47] + z[96]*q3p + z[97]*
  q2p - 2*w*z[32]*z[49];
  z[236] = z[38]*z[102] + 2*w*z[32]*z[48] + 2*w*z[34]*z[47] + z[100]*q3p - 
  0.5*z[36]*z[103] - 0.5*z[101]*q2p;
  z[237] = Izz*z[33]*z[231] + z[32]*(Ixx*z[232]+Ixy*z[233]) + m*(z[47]*z[234]+
  z[48]*z[235]+z[49]*z[236]) - z[34]*(Ixy*z[232]+Iyy*z[233]);

  no_cb[0] = q4 + z[26]*z[41] + z[39]*(z[9]+z[25]) - z[43]*(l+z[10]-z[27]);
  no_cb[1] = q5 + z[40]*(z[9]+z[25]) - z[26]*z[42] - z[44]*(l+z[10]-z[27]);
  no_cb[2] = -z[4]*z[26] - z[13]*(z[9]+z[25]) - z[14]*(l+z[10]-z[27]);
  H[0] = w*(Iyy*pow(z[34],2)+Izz*pow(z[33],2)+z[32]*(Ixx*z[32]-2*Ixy*z[34]));
  H[1] = -w*(Izz*z[4]*z[32]*z[33]-Izz*z[3]*z[6]*z[33]*z[34]-z[4]*z[33]*(Ixx*
  z[32]-Ixy*z[34])-z[3]*z[5]*z[34]*(Ixx*z[32]-Ixy*z[34])-(z[13]*z[33]+z[14]*
  z[32])*(Ixy*z[32]-Iyy*z[34]));
  H[2] = w*(Izz*z[14]*z[33]+z[4]*(Ixy*z[32]-Iyy*z[34])-z[13]*(Ixx*z[32]-Ixy*
  z[34]));
  p[0] = m*w*(z[32]*z[47]+z[33]*z[49]-z[34]*z[48]);
  p[1] = -m*w*(z[4]*z[32]*z[49]-z[4]*z[33]*z[47]-z[3]*z[5]*z[34]*z[47]-z[3]*
  z[6]*z[34]*z[49]-z[48]*(z[13]*z[33]+z[14]*z[32]));
  p[2] = -m*w*(z[13]*z[47]-z[4]*z[48]-z[14]*z[49]);
  df[0] = 0;
  df[1] = 0.5*w*z[120];
  df[2] = -0.5*w*z[129];
  df[3] = 0;
  df[4] = 0;
  df[5] = z[35];
  df[6] = 0;
  df[7] = w*z[130];
  df[8] = w*z[131];
  df[9] = 0;
  df[10] = 0;
  df[11] = z[36];
  df[12] = 0;
  df[13] = w*z[133];
  df[14] = w*z[135];
  df[15] = 0;
  df[16] = 0;
  df[17] = z[38];
  df[18] = -z[138];
  df[19] = -z[139];
  df[20] = -z[142];
  df[21] = 0;
  df[22] = 0;
  df[23] = -z[143];
  df[24] = -z[45];
  df[25] = -z[144];
  df[26] = -z[147];
  df[27] = 0;
  df[28] = 0;
  df[29] = -z[148];
  df[30] = 0;
  df[31] = z[188];
  df[32] = z[230];
  df[33] = 0;
  df[34] = 0;
  df[35] = z[237]/z[107];
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
