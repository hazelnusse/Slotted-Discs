autoz on
newtonian n
bodies S
frames a, b, db

constants m, ra, rb, k, l, Ix, Iy, Iz, g, alpha

variables q{3}', w1, w2
motionvariables u'
points ca, cb, mc

inertia S, Ix, Iy, Iz, 0, 0, 0
mass S=m

zee_not = [w1, w2, u, u', q1', q2', q3']
simprot(n, a, 3, q1)
simprot(a, b, 1, q2)
simprot(b, s, 2, q3)
simprot(s, db, 1, alpha)

p_ao_ca> = ra*unitvec(a3> - dot(a3>, s2>)*s2>)
p_ao_bo> = -l*s3>
p_bo_cb> = rb*unitvec(a3> - dot(a3>, db2>)*db2>)
p_ao_so> = -k*s3>

hc = dot(a3>, p_ca_cb>)

w_s_n> = w1*s1> + w2*s2> + u*s3>
kindiffs(n, s, BODY312, q1, q2, q3)
w_a_n> = q3'*a3>
vson1> = cross(w_s_n>, p_ca_so>)
vson2> = cross(w_s_n>, p_cb_so>)

nh[1] = dot(a1>, vson1> - vson2>)
nh[2] = dot(a2>, vson1> - vson2>)

solve(nh, [w1, w2])

q1' := replace(rhs(q1'), w1=rhs(w1), w2=rhs(w2))
q2' := replace(rhs(q2'), w1=rhs(w1), w2=rhs(w2))
q3' := replace(rhs(q3'), w1=rhs(w1), w2=rhs(w2))

w_s_n> := replace(w_s_n>, w1=rhs(w1), w2=rhs(w2))
v_so_n> = replace(vson1>, w1=rhs(w1), w2=rhs(w2))

zee_not := [u']

ke = m*mag(v_so_n>)/2
pe = -m*g*dot(p_ca_so>, a3>)

alf_s_n> = dt(w_s_n>, n)
a_so_n> = dt(v_so_n>, n)

fr_1 = dot(g*m*a3>, coef(v_so_n>, u))
fr_star_1 = -dot(m*a_so_n>, coef(v_so_n>, u)) - &
            dot(dot(alf_s_n>, I_S_SO>>) + cross(w_s_n>, dot(I_S_SO>>, w_s_n>)),  coef(w_s_n>, u))

solve(rhs(fr_1) + rhs(fr_star_1), u')

A[1] = d(rhs(q2'), q2)
A[2] = d(rhs(q2'), u)
A[3] = d(rhs(u'), q2)
A[4] = d(rhs(u'), u)

unitsystem  kg,m,s
input ra = 1.0 m, rb = 1.0 m, k = 0.5 m, l = 1.0 m
input alpha = pi/2 rad, m = 1.0 kg, g = 9.81 m/s^2
input Ix = 0.0 kg*m^2, Iy = 0.0 kg*m^2, Iz = 0.0 kg*m^2
input q1 = 0.0 rad, q2 = 0.0 rad, q3 = 0.0 %asin(rb/(ra+l))
input u = 0.0 rad/s

output q1 rad, q2 rad, q3 rad, q1' rad/s, q2' rad/s, q3' rad/s, u rad/s, ke kg*m^2/s/s, pe kg*m^2/s/s, hc m, A

code dynamics() slotted_discs_al.c


