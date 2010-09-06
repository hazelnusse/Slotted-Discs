autoz on
newtonian n
bodies da, db, s
frames a, b, c

constants m, ra, rb, k, l, g, alpha
constants Ixx, Iyy, Izz, Iyz

% Camera angles and distance from point in ground plane towards which it is
% pointed.
constants theta, phi, spin          % Camera inclination and azimuth
constants camx, camy, camz          % Camera target position
constants dcam                      % Distance from eye to camera target

points ct                           % Target point of camera

% Azimuth, Inclination, Camera spin
% Spin should be pi/2 to make c2> in a vertical plane
dircos(n, c, BODY323, phi, pi/2-theta, spin)

% Locating the inertial origin relative to the camera origin
p_co_ct> = -dcam*c3>
p_no_ct> = camx*n1> + camy*n2> + camz*n3>

variables q{5}'
variables w1', w2, w3
points ca, cb

mass s=m
inertia s, Ixx, Iyy, Izz, 0, Iyz, 0

simprot(n, a, 3, q1)
simprot(a, b, 2, q2)
simprot(b, da, 3, q3)
simprot(da, db, 1, alpha)
simprot(da, s, 1, 0)

% Disc A contact to Disc A center
p_ca_dao> = express(ra*unitvec(a3> - dot(a3>, da3>)*da3>), s)
% Disc A center to Disc B center
p_dao_dbo> = -l*s1>
% Disc A center to System Mass center
p_dao_so> = -k*s1>
% Disc B center to Disc B contact 
p_dbo_cb> = express(-rb*unitvec(a3> - dot(a3>, db3>)*db3>), s)
% Inertial origin to system mass center
p_no_so> = q4*n1> + q5*n2> + dot(p_ca_so>, a3>)*n3>

zee_not = [w1, w2, w3, q1', q2', q3']

% Angular velocity of system
w_s_n> = w1*s1> + w2*s2> + w3*s3>
% Velocity of disc B contact
vcbn> = cross(w_s_n>, p_ca_cb>)
nh[1] = dot(a1>, vcbn>)
nh[2] = dot(a2>, vcbn>)
solve(nh, [w2, w3])

kindiffs(N, S, BODY323, q1, q2, q3)
q1' := replace(q1', w2=rhs(w2), w3=rhs(w3))
q2' := replace(q2', w2=rhs(w2), w3=rhs(w3))
solve(dt(dot(p_ca_cb>, a3>)), q3')
q3' := replace(q3', q2'=rhs(q2'))

w_s_n> := replace(w_s_n>, w2=rhs(w2), w3=rhs(w3))
v_so_n> := cross(w_s_n>, p_ca_so>)

zee_not := [w1']
% Translation rates of system mass center
q4' = dot(v_so_n>, n1>)
q5' = dot(v_so_n>, n2>)

alf_s_n> = dt(w_s_n>, s)
a_so_n> = dt(v_so_n>, s) + cross(w_s_n>, v_so_n>)

fr_1 = dot(-g*m*a3>, coef(v_so_n>, w1))
fr_star_1 = -dot(m*a_so_n>, coef(v_so_n>, w1)) &
            -dot(dot(alf_s_n>, I_S_SO>>) + cross(w_s_n>, dot(I_S_SO>>, w_s_n>)), coef(w_s_n>, w1))
solve(rhs(fr_1) + rhs(fr_star_1), w1')

% Extraneous outputs
ke = m*dot(v_so_n>, v_so_n>)/2.0 + dot(w_s_n>, dot(I_S_SO>>, w_s_n>))/2.0
pe = m*g*dot(p_ca_so>, a3>)

T_da[1] = C_DA[1,1]
T_da[2] = C_DA[2,1]
T_da[3] = C_DA[3,1]
T_da[4] = 0
T_da[5] = C_DA[1,2]
T_da[6] = C_DA[2,2]
T_da[7] = C_DA[3,2]
T_da[8] = 0
T_da[9] = C_DA[1,3]
T_da[10] = C_DA[2,3]
T_da[11] = C_DA[3,3]
T_da[12] = 0
T_da[13] = dot(p_co_dao>, c1>)
T_da[14] = dot(p_co_dao>, c2>)
T_da[15] = dot(p_co_dao>, c3>)
T_da[16] = 1

T_db[1] = C_db[1,1]
T_db[2] = C_db[2,1]
T_db[3] = C_db[3,1]
T_db[4] = 0
T_db[5] = C_db[1,2]
T_db[6] = C_db[2,2]
T_db[7] = C_db[3,2]
T_db[8] = 0
T_db[9] = C_db[1,3]
T_db[10] = C_db[2,3]
T_db[11] = C_db[3,3]
T_db[12] = 0
T_db[13] = dot(p_co_dbo>, c1>)
T_db[14] = dot(p_co_dbo>, c2>)
T_db[15] = dot(p_co_dbo>, c3>)
T_db[16] = 1

% Contact point transformation matrices
T_ca[1] = dot(c1>, a1>)
T_ca[2] = dot(c2>, a1>)
T_ca[3] = dot(c3>, a1>)
T_ca[4] = 0
T_ca[5] = dot(c1>, a2>)
T_ca[6] = dot(c2>, a2>)
T_ca[7] = dot(c3>, a2>)
T_ca[8] = 0
T_ca[9] = dot(c1>, a3>)
T_ca[10] = dot(c2>, a3>)
T_ca[11] = dot(c3>, a3>)
T_ca[12] = 0
T_ca[13] = dot(p_co_ca>, c1>)
T_ca[14] = dot(p_co_ca>, c2>)
T_ca[15] = dot(p_co_ca>, c3>)
T_ca[16] = 1

% System mass center transformation matrix
T_so[1] = dot(c1>, s1>)
T_so[2] = dot(c2>, s1>)
T_so[3] = dot(c3>, s1>)
T_so[4] = 0
T_so[5] = dot(c1>, s2>)
T_so[6] = dot(c2>, s2>)
T_so[7] = dot(c3>, s2>)
T_so[8] = 0
T_so[9] = dot(c1>, s3>)
T_so[10] = dot(c2>, s3>)
T_so[11] = dot(c3>, s3>)
T_so[12] = 0
T_so[13] = dot(p_co_so>, c1>)
T_so[14] = dot(p_co_so>, c2>)
T_so[15] = dot(p_co_so>, c3>)
T_so[16] = 1

cb_tang> = unitvec(cross(p_dbo_cb>, a3>))
cb_perp> = cross(a3>, cb_tang>)

T_cb[1] = dot(c1>, cb_tang>)
T_cb[2] = dot(c2>, cb_tang>)
T_cb[3] = dot(c3>, cb_tang>)
T_cb[4] = 0
T_cb[5] = dot(c1>, cb_perp>)
T_cb[6] = dot(c2>, cb_perp>)
T_cb[7] = dot(c3>, cb_perp>)
T_cb[8] = 0
T_cb[9] = dot(c1>, a3>)
T_cb[10] = dot(c2>, a3>)
T_cb[11] = dot(c3>, a3>)
T_cb[12] = 0
T_cb[13] = dot(p_co_cb>, c1>)
T_cb[14] = dot(p_co_cb>, c2>)
T_cb[15] = dot(p_co_cb>, c3>)
T_cb[16] = 1

A[1] = d(rhs(q1'), q1)
A[2] = d(rhs(q1'), q2)
A[3] = d(rhs(q1'), q3)
A[4] = d(rhs(q1'), q4)
A[5] = d(rhs(q1'), q5)
A[6] = d(rhs(q1'), w1)

A[7] = d(rhs(q2'), q1)
A[8] = d(rhs(q2'), q2)
A[9] = d(rhs(q2'), q3)
A[10] = d(rhs(q2'), q4)
A[11] = d(rhs(q2'), q5)
A[12] = d(rhs(q2'), w1)

A[13] = d(rhs(q3'), q1)
A[14] = d(rhs(q3'), q2)
A[15] = d(rhs(q3'), q3)
A[16] = d(rhs(q3'), q4)
A[17] = d(rhs(q3'), q5)
A[18] = d(rhs(q3'), w1)

A[19] = d(rhs(q4'), q1)
A[20] = d(rhs(q4'), q2)
A[21] = d(rhs(q4'), q3)
A[22] = d(rhs(q4'), q4)
A[23] = d(rhs(q4'), q5)
A[24] = d(rhs(q4'), w1)

A[25] = d(rhs(q5'), q1)
A[26] = d(rhs(q5'), q2)
A[27] = d(rhs(q5'), q3)
A[28] = d(rhs(q5'), q4)
A[29] = d(rhs(q5'), q5)
A[30] = d(rhs(q5'), w1)

A[31] = d(rhs(q5'), q1)
A[32] = d(rhs(q5'), q2)
A[33] = d(rhs(q5'), q3)
A[34] = d(rhs(q5'), q4)
A[35] = d(rhs(q5'), q5)
A[36] = d(rhs(q5'), w1)

% Constraints, should be zero throughout numerical integration
con[1] = dot(p_ca_cb>, a3>)
con[2] = nh[1]
con[3] = nh[2]

% Linear momentum of system
p[1] = dot(m*v_so_n>, s1>)
p[2] = dot(m*v_so_n>, s2>)
p[3] = dot(m*v_so_n>, s3>)

% Angular momentum of system about mass center
Hvec> = dot(I_S_SO>>, w_s_n>)
H[1] = dot(Hvec>, s1>)
H[2] = dot(Hvec>, s2>)
H[3] = dot(Hvec>, s3>)


unitsystem  kg,m,s
input ra = 0.1 m, rb = 0.1 m, k = .05, l = .1 m
input alpha = pi/2 rad, m = 0.2 kg, g = 9.81 m/s^2
input Ixx = 1.0, Iyy = 1.0, Izz = 1.0, Iyz = 1.0
input q1 = 0.0 rad, q2 = 0.0 rad, q3 = 0.0, q4 = 0.0 m, q5 = 0.0 m
input w1 = 0.0 rad/s
input theta = pi/4.0, phi = pi/4.0, spin = pi/2.0
input camx = 0.0, camy = 0.0, camz = 0.0
input dcam = 1.0

output q1 rad, q2 rad, q3 rad, q4 m, q5 m, q1' rad/s, q2' rad/s, q3' rad/s
output q4' m/s, q5' m/s, w1 rad/s, w2 rad/s, w3 rad/s, ke kg*m^2/s/s, pe kg*m^2/s/s
encode A, T_da, T_db, T_ca, T_cb, T_so, con, p, H

code dynamics() slotted_discs_al.c
