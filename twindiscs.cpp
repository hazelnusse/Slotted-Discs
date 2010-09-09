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

SlottedDiscs::SlottedDiscs()
{
  // Set parameters
  ma = 2.0;
  mb = 2.0;
  ra = .1;
  rb = .1;
  l = sqrt(2.0)*ra;
  g = 9.81;
  alpha = M_PI/2.0;

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
  cout.precision(9);
  cout.width(18);
  cout << pe;
  cout.precision(9);
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
  cout << "ma = " << ma << endl << "mb = " << mb << endl;
  cout << "Ia = " << Ia << endl << "Ib = " << Ib << endl;
  cout << "Ja = " << Ja << endl << "Jb = " << Jb << endl;
  cout << "alpha = " << alpha << endl << "l = " << l << endl;
  cout << "g = " << g << endl;
} // printParameters()

void SlottedDiscs::evalConstants(void)
{
  Ib = 0.25*mb*pow(rb,2);
  z[7] = cos(alpha);
  Ia = 0.25*ma*pow(ra,2);
  Ja = 0.5*ma*pow(ra,2);
  Jb = 0.5*mb*pow(rb,2);
  z[8] = sin(alpha);
  z[157] = -Ia - Ib;
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
  z[71] = z[70]*(z[17]*(z[1]*z[6]+z[2]*z[10])+z[19]*(z[1]*z[5]-z[2]*z[9]));
  q4p = -w3*z[71];
  z[72] = z[70]*(z[19]*(z[1]*z[9]+z[2]*z[5])-z[17]*(z[1]*z[10]-z[2]*z[6]));
  q5p = -w3*z[72];
  z[41] = (z[33]*z[38]+z[35]*z[37])/z[39];
  z[73] = -z[7]*z[40] - z[8]*z[41];
  z[87] = z[10]*z[17] + z[12]*z[47] + z[10]*(z[7]*z[30]+z[8]*z[31]) - z[11]*(
  z[50]-z[65]) - z[9]*(l+z[19]-z[32]) - z[12]*(z[7]*z[59]+z[8]*z[63]);
  z[95] = z[12]*z[49] + z[3]*(l+z[19]-z[32]) + z[10]*(z[7]*z[31]-z[8]*z[30]) - 
  z[10]*z[18] - z[4]*(z[50]-z[65]) - z[12]*(z[7]*z[63]-z[8]*z[59]);
  z[84] = z[6]*z[47] + z[5]*(z[50]-z[65]) - z[6]*(z[7]*z[59]+z[8]*z[63]);
  z[94] = z[6]*(z[7]*z[63]-z[49]-z[8]*z[59]);
  z[97] = -z[33]*z[87] - z[34]*z[95] - z[35]*z[84] - z[36]*z[94];
  z[91] = z[5]*z[49] + z[21]*z[59] - z[20]*z[63];
  z[85] = z[4]*z[47] + z[9]*z[18] - z[3]*z[17] - z[11]*z[49] - z[22]*z[31] - 
  z[23]*z[30] - z[24]*z[63] - z[25]*z[59];
  z[99] = (z[97]*(z[34]*z[38]-z[36]*z[37])+z[39]*(z[36]*z[91]+z[37]*z[87]-
  z[34]*z[85]-z[38]*z[84]))/pow(z[39],2);
  z[103] = (z[97]*(z[33]*z[38]+z[35]*z[37])-z[39]*(z[33]*z[85]+z[35]*z[91]+
  z[37]*z[95]+z[38]*z[94]))/pow(z[39],2);
  z[111] = z[7]*z[99] + z[8]*z[103];
  z[113] = w3*z[111];
  z[88] = z[11]*z[66] + z[12]*(l-z[32]) + z[11]*(z[7]*z[30]+z[8]*z[31]) - 
  z[12]*(z[7]*z[60]+z[8]*z[64]);
  z[96] = z[11]*(z[7]*z[31]-z[8]*z[30]) - z[11]*z[18] - z[4]*(z[17]-z[66]) - 
  z[12]*(z[7]*z[64]-z[8]*z[60]);
  z[83] = z[6]*(l-z[32]) - z[5]*z[66] - z[5]*(z[7]*z[30]+z[8]*z[31]) - z[6]*(
  z[7]*z[60]+z[8]*z[64]);
  z[93] = z[6]*(z[7]*z[64]-z[8]*z[60]) + z[5]*(z[7]*z[31]-z[18]-z[8]*z[30]);
  z[98] = -z[33]*z[88] - z[34]*z[96] - z[35]*z[83] - z[36]*z[93];
  z[90] = z[7]*z[6];
  z[89] = z[8]*z[6];
  z[92] = z[21]*z[60] + z[31]*z[90] - z[6]*z[18] - z[20]*z[64] - z[30]*z[89];
  z[86] = z[31]*z[51] - z[4]*z[19] - z[12]*z[18] - z[24]*z[64] - z[25]*z[60] - 
  z[30]*z[52];
  z[100] = (z[98]*(z[34]*z[38]-z[36]*z[37])+z[39]*(z[36]*z[92]+z[37]*z[88]-
  z[34]*z[86]-z[38]*z[83]))/pow(z[39],2);
  z[104] = (z[98]*(z[33]*z[38]+z[35]*z[37])-z[39]*(z[33]*z[86]+z[35]*z[92]+
  z[37]*z[96]+z[38]*z[93]))/pow(z[39],2);
  z[112] = z[7]*z[100] + z[8]*z[104];
  z[114] = w3*z[112];
  z[119] = z[113]*q2p + z[114]*q3p;
  z[74] = z[8]*z[40] - z[7]*z[41];
  z[115] = z[7]*z[103] - z[8]*z[99];
  z[117] = w3*z[115];
  z[116] = z[7]*z[104] - z[8]*z[100];
  z[118] = w3*z[116];
  z[120] = z[117]*q2p + z[118]*q3p;
  z[75] = z[19]*z[41] - z[18];
  z[124] = -z[49] - z[19]*z[103] - z[41]*z[50];
  z[126] = w3*z[124];
  z[125] = -z[17]*z[41] - z[19]*z[104];
  z[127] = w3*z[125];
  z[76] = z[17] - z[19]*z[40];
  z[77] = z[17]*z[41] - z[18]*z[40];
  z[136] = z[126]*q2p + z[127]*q3p - pow(w3,2)*z[76] - pow(w3,2)*z[41]*z[77];
  z[128] = z[19] + z[17]*z[40] + z[19]*z[100];
  z[130] = w3*z[128];
  z[129] = z[19]*z[99] + z[40]*z[50] - z[47];
  z[131] = w3*z[129];
  z[137] = pow(w3,2)*z[75] + pow(w3,2)*z[40]*z[77] + z[130]*q3p + z[131]*q2p;
  z[132] = z[18]*z[99] - z[17]*z[103] - z[40]*z[49] - z[41]*z[47];
  z[134] = w3*z[132];
  z[133] = z[18]*z[100] + z[19]*z[41] - z[17]*z[104];
  z[135] = w3*z[133];
  z[138] = pow(w3,2)*z[41]*z[75] + z[134]*q2p + z[135]*q3p - pow(w3,2)*z[40]*
  z[76];
  z[78] = -z[31] - z[32]*z[74];
  z[80] = z[30]*z[74] + z[31]*z[73];
  z[142] = z[65]*z[74] - z[63] - z[32]*z[115];
  z[144] = w3*z[142];
  z[143] = z[66]*z[74] - z[64] - z[32]*z[116];
  z[145] = w3*z[143];
  z[79] = z[32]*z[73] - z[30];
  z[154] = pow(w3,2)*z[74]*z[80] + z[144]*q2p + z[145]*q3p - pow(w3,2)*z[79];
  z[146] = z[32]*z[111] - z[59] - z[65]*z[73];
  z[148] = w3*z[146];
  z[147] = z[32]*z[112] - z[60] - z[66]*z[73];
  z[149] = w3*z[147];
  z[155] = pow(w3,2)*z[78] + z[148]*q2p + z[149]*q3p - pow(w3,2)*z[73]*z[80];
  z[150] = z[30]*z[115] + z[31]*z[111] + z[59]*z[74] + z[63]*z[73];
  z[152] = w3*z[150];
  z[151] = z[30]*z[116] + z[31]*z[112] + z[60]*z[74] + z[64]*z[73];
  z[153] = w3*z[151];
  z[156] = pow(w3,2)*z[73]*z[79] + z[152]*q2p + z[153]*q3p - pow(w3,2)*z[74]*
  z[78];
  z[101] = w3*z[99];
  z[102] = w3*z[100];
  z[107] = z[101]*q2p + z[102]*q3p;
  z[105] = w3*z[103];
  z[106] = w3*z[104];
  z[108] = z[105]*q2p + z[106]*q3p;
  z[159] = Ib*z[73]*z[119] + Jb*z[74]*z[120] + ma*(z[75]*z[136]+z[76]*z[137]+
  z[77]*z[138]) + mb*(z[78]*z[154]+z[79]*z[155]+z[80]*z[156]) - g*(mb*(z[12]*
  z[80]+z[24]*z[78]+z[25]*z[79])-ma*(z[11]*z[75]-z[4]*z[76]-z[12]*z[77])) - 
  Ia*z[40]*z[107] - Ja*z[41]*z[108];
  z[158] = z[157] - Ia*pow(z[40],2) - Ib*pow(z[73],2) - Ja*pow(z[41],2) - Jb*
  pow(z[74],2) - ma*(pow(z[75],2)+pow(z[76],2)+pow(z[77],2)) - mb*(pow(z[78],
  2)+pow(z[79],2)+pow(z[80],2));
  z[160] = z[159]/z[158];
  w3p = z[160];
} // eoms()

void SlottedDiscs::computeOutputs(void)
{
  w1 = -w3*z[40];
  w2 = -w3*z[41];
  ke = 0.5*pow(w3,2)*(Ia+Ib+Ia*pow(z[40],2)+Ib*pow(z[73],2)+Ja*pow(z[41],2)+
  Jb*pow(z[74],2)+ma*(pow(z[75],2)+pow(z[76],2)+pow(z[77],2))+mb*(pow(z[78],2)+
  pow(z[79],2)+pow(z[80],2)));
  pe = -g*(ma*(z[4]*z[18]-z[11]*z[17]-z[12]*z[19])+mb*(z[4]*z[18]-z[11]*z[17]-
  z[12]*(l+z[19])));
  te = ke + pe;
  z[161] = z[1]*z[5] - z[2]*z[9];
  z[162] = z[1]*z[9] + z[2]*z[5];
  z[163] = z[2]*z[3];
  z[164] = z[1]*z[3];
  z[165] = z[1]*z[6] + z[2]*z[10];
  z[166] = z[2]*z[6] - z[1]*z[10];
  z[167] = z[7]*z[161] - z[8]*z[163];
  z[168] = z[7]*z[162] + z[8]*z[164];
  z[169] = -z[7]*z[163] - z[8]*z[161];
  z[170] = z[7]*z[164] - z[8]*z[162];

  con[0] = w3*z[20]*z[31] + z[6]*(w1*z[18]-w2*z[17]) + z[5]*(w2*z[32]-w3*
  z[18]-w2*(l+z[19])) + z[6]*(z[8]*w1*z[30]-z[7]*w1*z[31]-z[7]*w2*z[30]-z[8]*
  w2*z[31]) - w3*z[21]*z[30];
  con[1] = w3*z[24]*z[31] + w3*z[25]*z[30] + z[12]*(w1*z[18]-w2*z[17]) + 
  z[12]*(z[8]*w1*z[30]-z[7]*w1*z[31]-z[7]*w2*z[30]-z[8]*w2*z[31]) - z[4]*(w1*
  z[32]-w3*z[17]-w1*(l+z[19])) - z[11]*(w2*z[32]-w3*z[18]-w2*(l+z[19]));
  con[2] = z[4]*z[18] + z[12]*z[32] + z[24]*z[30] - z[11]*z[17] - z[25]*z[31] - 
  z[12]*(l+z[19]);
  no_cb[0] = q4 + z[17]*z[161] + z[30]*z[167] + z[32]*z[165] - z[18]*z[163] - 
  z[31]*z[169] - z[165]*(l+z[19]);
  no_cb[1] = q5 + z[17]*z[162] + z[18]*z[164] + z[30]*z[168] + z[32]*z[166] - 
  z[31]*z[170] - z[166]*(l+z[19]);
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
  ma = p->ma;
  ra = p->ra;
  mb = p->mb;
  mb = p->mb;
  l = p->l;
  g = p->g;
  alpha = p->alpha;
  // Constants
  Ib = 0.25*mb*pow(rb,2);
  z[7] = cos(alpha);
  Ia = 0.25*ma*pow(ra,2);
  Ja = 0.5*ma*pow(ra,2);
  Jb = 0.5*mb*pow(rb,2);
  z[8] = sin(alpha);
  z[157] = -Ia - Ib;
} // setParameters()
