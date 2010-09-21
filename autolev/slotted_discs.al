newtonian n
bodies s
frames da, db
frames a, b

constants m, ra, rb, l, g, alpha
constants k

variables q{5}'
variables w'
points ca, cb

inertia s(da), Ixx, Iyy, Izz, Ixy, 0, 0

autoz on
simprot(n, a, 3, q1)
simprot(a, b, 1, q2)
simprot(b, da, 2, q3)
simprot(da, db, 3, alpha)

p_no_ca> = q4*n1> + q5*n2>
p_ca_dao> = express(-ra*b3>, da)
p_dao_dbo> = -l*da3>
p_dao_so> = -k*da3>
p_dbo_cb> = express(rb*unitvec(a3> - dot(a3>, db2>)*db2>), da)

zee_not = [w]

cl1> = unitvec(p_ca_cb>)
cl2> = cross(a3>, cl1>)

w_da_n> = w*cl1>
w_db_da> = 0>

kindiffs(N, DA, BODY312, q1, q2, q3)

zero> = q4'*n1> + q5'*n2> - cross(q3'*b2>, p_ca_dao>)
wr[1] = dot(zero>, n1>)
wr[2] = dot(zero>, n2>)
solve(wr, [q4', q5'])

v_so_n> = cross(w_da_n>, p_ca_so>)
zee_not := [w']

alf_da_n> = dt(w_da_n>, da)
a_so_n> = dt(v_so_n>, da) + cross(w_da_n>, v_so_n>)

fr_1 = dot(g*m*a3>, coef(v_so_n>, w))
fr_star_1 = -dot(m*a_so_n>, coef(v_so_n>, w)) &
            -dot(dot(alf_da_n>, I_S_SO>>) + cross(w_da_n>, dot(I_S_SO>>, w_da_n>)), coef(w_da_n>, w))

solve(rhs(fr_1) + rhs(fr_star_1), w')

% Extraneous outputs
ke = m*dot(v_so_n>, v_so_n>)/2.0 + dot(w_da_n>, dot(I_S_SO>>, w_da_n>))/2.0
pe = -m*g*dot(p_ca_so>, a3>)
te = ke + pe

no_cb[1] = dot(p_no_cb>, n1>)
no_cb[2] = dot(p_no_cb>, n2>)
no_cb[3] = dot(p_no_cb>, a3>)

w1 = dot(w_da_n>, da1>)
w2 = dot(w_da_n>, da2>)
w3 = dot(w_da_n>, da3>)

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
input ra = 0.1 m, rb = 0.1 m, l = .1 m, k = 0.0 m
input alpha = pi/2 rad, m = 0.1 kg, g = 9.81 m/s^2
input Ixx = 0.0, Iyy = 0.0, Izz = 0.0, Ixy = 0.0
input q1 = 0.0 rad, q2 = 0.0 rad, q3 = 0.0, q4 = 0.0 m, q5 = 0.0 m
input w = 0.0 rad/s

output q1 rad, q2 rad, q3 rad, q4 m, q5 m, q1' rad/s, q2' rad/s, q3' rad/s
output q4' m/s, q5' m/s, w1 rad/s, w2 rad/s, w3 rad/s, ke kg*m^2/s/s, pe kg*m^2/s/s, te kg*m^2/s/s
encode no_cb, H, p, df

code dynamics() slotted_discs_al.c
