newtonian n
bodies da, db
frames a, b

constants ma, mb, ra, rb, l, g, alpha

variables q{5}'
variables w1, w2, w3'
points ca, cb

Ia = ma*ra*ra/4.0
Ja = ma*ra*ra/2.0

Ib = mb*rb*rb/4.0
Jb = mb*rb*rb/2.0

inertia da, Ia, Ja, Ia, 0, 0, 0
inertia db, Ib, Jb, Ib, 0, 0, 0

autoz on
simprot(n, a, 3, q1)
simprot(a, b, 1, q2)
simprot(b, da, 2, q3)
simprot(da, db, 3, alpha)

p_no_ca> = q4*n1> + q5*n2>
p_ca_dao> = express(-ra*unitvec(a3> - dot(a3>, da2>)*da2>), da)
p_dao_dbo> = -l*da3>
p_dbo_cb> = express(rb*unitvec(a3> - dot(a3>, db2>)*db2>), db)

autoz off
hc = explicit(dot(p_ca_cb>, a3>))
hc_df = explicit(d(hc, q3))
autoz on

zee_not = [w1, w2, w3, q1', q2', q3']
w_da_n> = w1*da1> + w2*da2> + w3*da3>
w_db_n> = express(w_da_n>, db)
% Velocity of cb
vcbn> = cross(w_da_n>, p_ca_cb>)
nh[1] = dot(a1>, vcbn>)
nh[2] = dot(a3>, vcbn>)
solve(nh, [w1, w2])

q1' = zee((-rhs(w1)*sin(q3) + w3*cos(q3))/cos(q2))
q2' = zee(rhs(w1)*cos(q3) + w3*sin(q3))
solve(dt(dot(p_ca_cb>, a3>)), q3')
q3' := zee(replace(q3', q2'=rhs(q2')))
% q3' = (rhs(w1)*sin(q3) - w3*cos(q3))*tan(q2) + rhs(w2)

autoz off
vcon1> = q4'*n1> + q5'*n2> + cross(q1'*a3> + q2'*a1>, p_ca_dao>)
vcon2> = cross(q1'*a3> + q2'*a1> + q3'*b2>, p_ca_dao>)
wr[1] = dot(vcon1> - vcon2>, n1>)
wr[2] = dot(vcon1> - vcon2>, n2>)
solve(wr, [q4', q5'])
q4' := zee(replace(q4', q3'=rhs(q3')))
q5' := zee(replace(q5', q3'=rhs(q3')))
autoz on

w_da_n> := replace(w_da_n>, w1=rhs(w1), w2=rhs(w2))
w_db_n> := replace(w_db_n>, w1=rhs(w1), w2=rhs(w2))
v_dao_n> = cross(w_da_n>, p_ca_dao>)
v_dbo_n> = cross(w_db_n>, p_cb_dbo>)

zee_not := [w3']

alf_da_n> = dt(w_da_n>, da)
alf_db_n> = dt(w_db_n>, db)
a_dao_n> = dt(v_dao_n>, da) + cross(w_da_n>, v_dao_n>)
a_dbo_n> = dt(v_dbo_n>, db) + cross(w_db_n>, v_dbo_n>)

fr_1 = dot(g*ma*a3>, coef(v_dao_n>, w3)) + dot(g*mb*a3>, coef(v_dbo_n>, w3))
fr_star_1 = -dot(ma*a_dao_n>, coef(v_dao_n>, w3)) &
            -dot(mb*a_dbo_n>, coef(v_dbo_n>, w3)) &
            -dot(dot(alf_da_n>, I_DA_DAO>>) + cross(w_da_n>, dot(I_DA_DAO>>, w_da_n>)), coef(w_da_n>, w3)) &
            -dot(dot(alf_db_n>, I_DB_DBO>>) + cross(w_db_n>, dot(I_DB_DBO>>, w_db_n>)), coef(w_db_n>, w3))
pause
solve(rhs(fr_1) + rhs(fr_star_1), w3')


% Extraneous outputs
ke = ma*dot(v_dao_n>, v_dao_n>)/2.0 + &
     mb*dot(v_dbo_n>, v_dbo_n>)/2.0 + &
     dot(w_da_n>, dot(I_DA_DAO>>, w_da_n>))/2.0 + &
     dot(w_db_n>, dot(I_DB_DBO>>, w_db_n>))/2.0
pe = - ma*g*dot(p_ca_dao>, a3>) - mb*g*dot(p_ca_dbo>, a3>)

T_da[1] = N_DA[1,1]
T_da[2] = N_DA[2,1]
T_da[3] = N_DA[3,1]
T_da[4] = 0
T_da[5] = N_DA[1,2]
T_da[6] = N_DA[2,2]
T_da[7] = N_DA[3,2]
T_da[8] = 0
T_da[9] = N_DA[1,3]
T_da[10] = N_DA[2,3]
T_da[11] = N_DA[3,3]
T_da[12] = 0
T_da[13] = dot(p_no_dao>, n1>)
T_da[14] = dot(p_no_dao>, n2>)
T_da[15] = dot(p_no_dao>, n3>)
T_da[16] = 1

T_db[1] = N_db[1,1]
T_db[2] = N_db[2,1]
T_db[3] = N_db[3,1]
T_db[4] = 0
T_db[5] = N_db[1,2]
T_db[6] = N_db[2,2]
T_db[7] = N_db[3,2]
T_db[8] = 0
T_db[9] = N_db[1,3]
T_db[10] = N_db[2,3]
T_db[11] = N_db[3,3]
T_db[12] = 0
T_db[13] = dot(p_no_dbo>, n1>)
T_db[14] = dot(p_no_dbo>, n2>)
T_db[15] = dot(p_no_dbo>, n3>)
T_db[16] = 1

no_cb[1] = dot(p_no_cb>, n1>)
no_cb[2] = dot(p_no_cb>, n2>)
no_cb[3] = dot(p_no_cb>, a3>)

states = [q1, q2, q3, q4, q5, w3]
A = d([rhs(q1'); rhs(q2'); rhs(q3'); rhs(q4'); rhs(q5'); rhs(w3')], states)

unitsystem  kg,m,s
input ra = 0.1 m, rb = 0.1 m, l = .1 m
input alpha = pi/2 rad, ma = 0.1 kg, mb = 0.1 kg, g = 9.81 m/s^2
input q1 = 0.0 rad, q2 = 0.0 rad, q3 = 0.0, q4 = 0.0 m, q5 = 0.0 m
input w3 = 0.0 rad/s

output q1 rad, q2 rad, q3 rad, q4 m, q5 m, q1' rad/s, q2' rad/s, q3' rad/s
output q4' m/s, q5' m/s, w1 rad/s, w2 rad/s, w3 rad/s, ke kg*m^2/s/s, pe kg*m^2/s/s
encode A, no_cb, T_da, T_db

code dynamics() slotted_discs_al.c
