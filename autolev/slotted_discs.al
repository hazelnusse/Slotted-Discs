autoz on
autorhs on

newtonian n
bodies s
frames da, db, dagl, dbgl
frames a, b, cam

constants m, ra, rb, l, g, alpha
constants k

variables q{5}'
variables w'
variables fx, fay, faz, fby, fbz
constants au1, au2, au3, au4, au5, au6

points ca, cb

% Camera variables and settings
points ct        % camera target
% phi: elevation from horizontal
% theta: azimuth
% d:    distance from camera targe along camera z
% ctx, cty, ctz:  x, y, z coorindates of camera target in inertial coordinates
constants phi, theta, d, ctx, cty, ctz

dircos(cam, n, body121, theta, phi, pi/2)
p_camo_ct> = -d*cam3>
p_no_ct> = ctx*n1> + cty*n2> + ctz*n3>

inertia s(da), Ixx, Iyy, Izz, Ixy, 0, 0
mass s=m

simprot(n, a, 3, q1)
simprot(a, b, 1, q2)
simprot(b, da, 2, q3)
simprot(da, db, 3, alpha)

% OpenGL functions draw most symmetric objects with axis of symmetry along
%z-axis, but since our axis of symmetry is the y axis, we need to do a rotation
%about the x-axis by -pi/2
simprot(da, dagl, 1, -pi/2)
simprot(db, dbgl, 1, -pi/2)

p_no_ca> = q4*n1> + q5*n2>
p_ca_dao> = express(-ra*b3>, da)
p_dao_dbo> = -l*da3>
p_dao_so> = -k*da3>
p_dbo_cb> = express(rb*unitvec(a3> - dot(a3>, db2>)*db2>), da)

zee_not = [w]

cl1> = unitvec(p_ca_cb>)
cl2> = cross(a3>, cl1>)

% Disc B contact
cb1> = cross(a3>, unitvec(p_dbo_cb>))
cb2> = cross(a3>, cb1>)

w_da_n> = w*cl1>
w_db_da> = 0>

w1 = dot(w_da_n>, da1>)
w2 = dot(w_da_n>, da2>)
w3 = dot(w_da_n>, da3>)

q1' = (-w1*zee(sin(q3)) + w3*zee(cos(q3)))/zee(cos(q2))
q2' = w1*zee(cos(q3)) + w3*zee(sin(q3))
autoz off
solve(dt(dot(p_ca_cb>, a3>)), q3')
zee(q3')
autoz on

zero> = q4'*n1> + q5'*n2> - cross(q3'*b2>, p_ca_dao>)
wr[1] = dot(zero>, n1>)
wr[2] = dot(zero>, n2>)
solve(wr, [q4', q5'])

zee_not := [w, fx, fay, faz, fby, fbz]
v_so_n> = cross(w*cl1>, p_ca_so>)

zee_not := [w', fx, fay, faz, fby, fbz]

alf_da_n> = dt(w_da_n>, da)
a_so_n> = dt(v_so_n>, da) + cross(w_da_n>, v_so_n>)

f_so> = g*m*a3>

R_star> = - m*a_so_n>
T_star> = - dot(alf_da_n>, I_S_SO>>) - cross(w_da_n>, dot(I_S_SO>>, w_da_n>))

% Single motion equation
solve(dot(f_so>, cross(cl1>, p_ca_so>)) + dot(cl1>, T_star>) + dot(cross(cl1>, p_ca_so>), R_star>), w')

% Constraint force determination
f_ca> = fay*cl2> + faz*a3>
f_cb> = fby*cl2> + fbz*a3>

% Moments about Disc A contact
t_ca> = cross(p_ca_cb>, f_cb>) + cross(p_ca_so>, f_so>)
H_ca> = dot(I_S_CA>>, w_da_n>)
H_ca_dot> = dt(H_ca>, da) + cross(w_da_n>, H_ca>)
H_ca_dot> := H_ca_dot>
NE1> = t_ca> - H_ca_dot>
con_eqns1 = [dot(NE1>, cl2>), dot(NE1>, a3>)]
solve(con_eqns1, [fby, fbz])

% Moments about Disc B contact
t_cb> = cross(p_cb_ca>, f_ca>) + cross(p_cb_so>, f_so>)
H_cb> = dot(I_S_CB>>, w_da_n>)
H_cb_dot> = dt(H_cb>, da) + cross(w_da_n>, H_cb>)
H_cb_dot> := H_cb_dot>
NE2> = t_cb> - H_cb_dot>
con_eqns2 = [dot(NE2>, cl2>), dot(NE2>, a3>)]
solve(con_eqns2, [fay, faz])

% Forces along contact line are indeterminate, best we can do is solve for
% resultant along contact line
fx = m*dot(a_so_n>, cl1>)

% Forward dynamics problem -- given the forces calculated above, we should be
% able to apply them to another rigid body with identical geometry and inertia
% and it should move identically.  This will be used as a check that the
% constraint forces were calculated correctly
%body s2
%frames da2, db2
%frames a2, b2
%points ca2, cb2
%inertia s2(da2), Ixx, Iyy, Izz, Ixy, 0, 0
%mass s2=m

%variables q2{5}'
%simprot(n, a2, 3, q21)
%simprot(a2, b2, 1, q22)
%simprot(b2, da2, 2, q23)
%simprot(da2, db2, 3, alpha)

%p_no_ca2> = q24*n1> + q25*n2>
%p_ca2_da2o> = express(-ra*b23>, da2)
%p_da2o_db2o> = -l*da23>
%p_da2o_s2o> = -k*da23>
%p_db2o_cb2> = express(rb*unitvec(a23> - dot(a23>, db2>)*db2>), da2)

% Extraneous outputs
ke = m*dot(v_so_n>, v_so_n>)/2.0 + dot(w_da_n>, dot(I_S_SO>>, w_da_n>))/2.0
pe = -m*g*dot(p_ca_so>, a3>)
te = ke + pe

no_cb[1] = dot(p_no_cb>, n1>)
no_cb[2] = dot(p_no_cb>, n2>)
no_cb[3] = dot(p_no_cb>, a3>)

no_so[1] = dot(p_no_so>, n1>)
no_so[2] = dot(p_no_so>, n2>)
no_so[3] = dot(p_no_so>, a3>)

a_so[1] = dot(a_so_n>, cl1>)
a_so[2] = dot(a_so_n>, cl2>)
a_so[3] = dot(a_so_n>, a3>)

% Angular momentum of system about system mass center
H_SYS_SO> = dot(I_S_SO>>, w_da_n>)
% Resolve H into components of the contact line coordinate system
H[1] = dot(H_SYS_SO>, cl1>)
H[2] = dot(H_SYS_SO>, cl2>)
H[3] = dot(H_SYS_SO>, a3>)

% Linear momentum of system mass center
p[1] = m*dot(v_so_n>, cl1>)
p[2] = m*dot(v_so_n>, cl2>)
p[3] = m*dot(v_so_n>, a3>)

% 4x4 OpenGL transformation matrices as column major arrays

% Disc A center
T_da[1] = dot(cam1>, da1>)
T_da[2] = dot(cam2>, da1>)
T_da[3] = dot(cam3>, da1>)
T_da[4] = 0
T_da[5] = dot(cam1>, da2>)
T_da[6] = dot(cam2>, da2>)
T_da[7] = dot(cam3>, da2>)
T_da[8] = 0
T_da[9] = dot(cam1>, da3>)
T_da[10] = dot(cam2>, da3>)
T_da[11] = dot(cam3>, da3>)
T_da[12] = 0
T_da[13] = dot(p_camo_dao>, cam1>)
T_da[14] = dot(p_camo_dao>, cam2>)
T_da[15] = dot(p_camo_dao>, cam3>)
T_da[16] = 1

% Disc B center
T_db[1] = dot(cam1>, db1>)
T_db[2] = dot(cam2>, db1>)
T_db[3] = dot(cam3>, db1>)
T_db[4] = 0
T_db[5] = dot(cam1>, db2>)
T_db[6] = dot(cam2>, db2>)
T_db[7] = dot(cam3>, db2>)
T_db[8] = 0
T_db[9] = dot(cam1>, db3>)
T_db[10] = dot(cam2>, db3>)
T_db[11] = dot(cam3>, db3>)
T_db[12] = 0
T_db[13] = dot(p_camo_dbo>, cam1>)
T_db[14] = dot(p_camo_dbo>, cam2>)
T_db[15] = dot(p_camo_dbo>, cam3>)
T_db[16] = 1

% Disc A center for OpenGL shapes
T_dagl[1] = dot(cam1>, dagl1>)
T_dagl[2] = dot(cam2>, dagl1>)
T_dagl[3] = dot(cam3>, dagl1>)
T_dagl[4] = 0
T_dagl[5] = dot(cam1>, dagl2>)
T_dagl[6] = dot(cam2>, dagl2>)
T_dagl[7] = dot(cam3>, dagl2>)
T_dagl[8] = 0
T_dagl[9] = dot(cam1>, dagl3>)
T_dagl[10] = dot(cam2>, dagl3>)
T_dagl[11] = dot(cam3>, dagl3>)
T_dagl[12] = 0
T_dagl[13] = dot(p_camo_dao>, cam1>)
T_dagl[14] = dot(p_camo_dao>, cam2>)
T_dagl[15] = dot(p_camo_dao>, cam3>)
T_dagl[16] = 1

% Disc B center for OpenGL shapes
T_dbgl[1] = dot(cam1>, dbgl1>)
T_dbgl[2] = dot(cam2>, dbgl1>)
T_dbgl[3] = dot(cam3>, dbgl1>)
T_dbgl[4] = 0
T_dbgl[5] = dot(cam1>, dbgl2>)
T_dbgl[6] = dot(cam2>, dbgl2>)
T_dbgl[7] = dot(cam3>, dbgl2>)
T_dbgl[8] = 0
T_dbgl[9] = dot(cam1>, dbgl3>)
T_dbgl[10] = dot(cam2>, dbgl3>)
T_dbgl[11] = dot(cam3>, dbgl3>)
T_dbgl[12] = 0
T_dbgl[13] = dot(p_camo_dbo>, cam1>)
T_dbgl[14] = dot(p_camo_dbo>, cam2>)
T_dbgl[15] = dot(p_camo_dbo>, cam3>)
T_dbgl[16] = 1

% Center of mass
T_so[1] = dot(cam1>, da1>)
T_so[2] = dot(cam2>, da1>)
T_so[3] = dot(cam3>, da1>)
T_so[4] = 0
T_so[5] = dot(cam1>, da2>)
T_so[6] = dot(cam2>, da2>)
T_so[7] = dot(cam3>, da2>)
T_so[8] = 0
T_so[9] = dot(cam1>, da3>)
T_so[10] = dot(cam2>, da3>)
T_so[11] = dot(cam3>, da3>)
T_so[12] = 0
T_so[13] = dot(p_camo_so>, cam1>)
T_so[14] = dot(p_camo_so>, cam2>)
T_so[15] = dot(p_camo_so>, cam3>)
T_so[16] = 1

% Disc A contact
T_ca[1] = dot(cam1>, a1>)
T_ca[2] = dot(cam2>, a1>)
T_ca[3] = dot(cam3>, a1>)
T_ca[4] = 0
T_ca[5] = dot(cam1>, a2>)
T_ca[6] = dot(cam2>, a2>)
T_ca[7] = dot(cam3>, a2>)
T_ca[8] = 0
T_ca[9] = dot(cam1>, a3>)
T_ca[10] = dot(cam2>, a3>)
T_ca[11] = dot(cam3>, a3>)
T_ca[12] = 0
T_ca[13] = dot(p_camo_ca>, cam1>)
T_ca[14] = dot(p_camo_ca>, cam2>)
T_ca[15] = dot(p_camo_ca>, cam3>)
T_ca[16] = 1

T_cb[1] = dot(cam1>, cb1>)
T_cb[2] = dot(cam2>, cb1>)
T_cb[3] = dot(cam3>, cb1>)
T_cb[4] = 0
T_cb[5] = dot(cam1>, cb2>)
T_cb[6] = dot(cam2>, cb2>)
T_cb[7] = dot(cam3>, cb2>)
T_cb[8] = 0
T_cb[9] = dot(cam1>, a3>)
T_cb[10] = dot(cam2>, a3>)
T_cb[11] = dot(cam3>, a3>)
T_cb[12] = 0
T_cb[13] = dot(p_camo_cb>, cam1>)
T_cb[14] = dot(p_camo_cb>, cam2>)
T_cb[15] = dot(p_camo_cb>, cam3>)
T_cb[16] = 1

% Jacobian of system of 6 ODE's, stored in a length 36 row-major array
df[1] = d(q1', q1)
df[2] = d(q1', q2)
df[3] = d(q1', q3)
df[4] = d(q1', q4)
df[5] = d(q1', q5)
df[6] = d(q1', w)

df[7] = d(q2', q1)
df[8] = d(q2', q2)
df[9] = d(q2', q3)
df[10] = d(q2', q4)
df[11] = d(q2', q5)
df[12] = d(q2', w)

df[13] = d(q3', q1)
df[14] = d(q3', q2)
df[15] = d(q3', q3)
df[16] = d(q3', q4)
df[17] = d(q3', q5)
df[18] = d(q3', w)

df[19] = d(q4', q1)
df[20] = d(q4', q2)
df[21] = d(q4', q3)
df[22] = d(q4', q4)
df[23] = d(q4', q5)
df[24] = d(q4', w)

df[25] = d(q5', q1)
df[26] = d(q5', q2)
df[27] = d(q5', q3)
df[28] = d(q5', q4)
df[29] = d(q5', q5)
df[30] = d(q5', w)

df[31] = d(w', q1)
df[32] = d(w', q2)
df[33] = d(w', q3)
df[34] = d(w', q4)
df[35] = d(w', q5)
df[36] = d(w', w)

unitsystem  kg,m,s
input theta = pi/4, phi = 0.0, ctx = 0.0, cty = 0.0, ctz = 0.0
input ra = 0.1 m, rb = 0.1 m, l = .1 m, k = 0.0 m
input alpha = pi/2 rad, m = 0.1 kg, g = 9.81 m/s^2
input Ixx = 0.0, Iyy = 0.0, Izz = 0.0, Ixy = 0.0
input q1 = 0.0 rad, q2 = 0.0 rad, q3 = 0.0, q4 = 0.0 m, q5 = 0.0 m
input w = 0.0 rad/s, d = 1.5 m

output q1 rad, q2 rad, q3 rad, q4 m, q5 m, q1' rad/s, q2' rad/s, q3' rad/s
output q4' m/s, q5' m/s, w1 rad/s, w2 rad/s, w3 rad/s, ke kg*m^2/s/s, pe kg*m^2/s/s, te kg*m^2/s/s
output fx, fay, faz, fby, fbz
encode no_cb, no_so, a_so, H, p, df, T_da, T_db, T_so, T_ca, T_cb, T_dagl, T_dbgl

code dynamics() slotted_discs_al.c
