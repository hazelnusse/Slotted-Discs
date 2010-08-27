autoz on
newtonian n
bodies S
frames a, b, db

constants m, ra, rb, k, l, Ix, Iy, Iz, g, alpha

variables q{5}'
specified w{2}'
motionvariables u'
points ca, cb

inertia S, Ix, Iy, Iz, 0, 0, 0
mass S=m

zee_not = [w1, w2, u, u', q1', q2', q3', q4', q5']
simprot(n, a, 3, q1)
simprot(a, b, 1, q2)
simprot(b, s, 2, q3)
simprot(s, db, 3, alpha)

p_no_ca> = q4*n1> + q5*n2>
p_ao_ca> = ra*unitvec(a3> - dot(a3>, s2>)*s2>)
p_ao_bo> = -l*s3>
p_bo_cb> = rb*unitvec(a3> - dot(a3>, db2>)*db2>)
p_ao_so> = -k*s3>

w_s_n> = w1*s1> + w2*s2> + u*s3>
q1' = (-w1*sin(q3) + u*cos(q3))/cos(q2)
solve(dt(dot(p_ca_cb>, a3>)), q2')
q3' = (w1*sin(q3) - u*cos(q3))*tan(q2) + w2
w_a_n> = q3'*a3>

autoz off
vcon1> = q4'*n1> + q5'*n2> + cross(q1'*a3> + q2'*a1>, p_ca_ao>)
vcon2> = cross(q1'*a3> + q2'*a1> + q3'*b2>, p_ca_ao>)
wr[1] = dot(vcon1> - vcon2>, n1>)
wr[2] = dot(vcon1> - vcon2>, n2>)
solve(wr, [q4', q5'])
autoz on

vson1> = cross(w_s_n>, p_ca_so>)
vson2> = cross(w_s_n>, p_cb_so>)

nh[1] = dot(a1>, vson1> - vson2>)
nh[2] = dot(a3>, vson1> - vson2>)

solve(nh, [w1, w2])

w_s_n> := replace(w_s_n>, w1=rhs(w1), w2=rhs(w2))
v_so_n> = cross(w_s_n>, p_ca_so>)

zee_not := [u']

alf_s_n> = dt(w_s_n>, n)
a_so_n> = dt(v_so_n>, n)

pause
fr_1 = dot(g*m*a3>, coef(v_so_n>, u))
fr_star_1 = -dot(m*a_so_n>, coef(v_so_n>, u)) - &
            dot(dot(alf_s_n>, I_S_SO>>) + cross(w_s_n>, dot(I_S_SO>>, w_s_n>)),  coef(w_s_n>, u))
solve(rhs(fr_1) + rhs(fr_star_1), u')
% Extraneous outputs
ke = m*mag(v_so_n>)/2
pe = -m*g*dot(p_ca_so>, a3>)
T_da = [N_S[1,1], N_S[2,1], N_S[3,1], 0, &
        N_S[1,2], N_S[2,2], N_S[3,2], 0, &
        N_S[1,3], N_S[2,3], N_S[3,3], 0, &
        dot(p_no_ao>, n1>), dot(p_no_ao>, n2>), dot(p_no_ao>, n3>), 1]

T_db = [N_db[1,1], N_db[2,1], N_db[3,1], 0, &
        N_db[1,2], N_db[2,2], N_db[3,2], 0, &
        N_db[1,3], N_db[2,3], N_db[3,3], 0, &
        dot(p_no_bo>, n1>), dot(p_no_bo>, n2>), dot(p_no_bo>, n3>), 1]

no_cb = [dot(p_no_cb>, n1>), dot(p_no_cb>, n2>), dot(p_no_cb>, a3>)]

A[1] = d(rhs(q2'), q2)
A[2] = d(rhs(q2'), u)
A[3] = d(rhs(u'), q2)
A[4] = d(rhs(u'), u)

unitsystem  kg,m,s
input ra = 1.0 m, rb = 1.0 m, k = 0.5 m, l = 1.0 m
input alpha = pi/2 rad, m = 1.0 kg, g = 9.81 m/s^2
input Ix = 0.0 kg*m^2, Iy = 0.0 kg*m^2, Iz = 0.0 kg*m^2
input q1 = 0.0 rad, q2 = 0.0 rad, q3 = 0.0, q4 = 0.0 m, q5 = 0.0 m %asin(rb/(ra+l))
input u = 0.0 rad/s

output q1 rad, q2 rad, q3 rad, q4 m, q5 m, q1' rad/s, q2' rad/s, q3' rad/s, q4' m/s, q5' m/s, u rad/s, ke kg*m^2/s/s, pe kg*m^2/s/s, A, no_cb

code dynamics() slotted_discs_al.c
