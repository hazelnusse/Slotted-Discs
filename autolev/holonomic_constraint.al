autoz off
newtonian n
frames a, b, da, db

constants m, ra, rb, l, g, alpha
constants k

variables q0, q1, q2

points ca, cb

simprot(n, a, 3, q0)
simprot(a, b, 1, q1)
simprot(b, da, 2, q2)
simprot(da, db, 3, alpha)

p_ca_dao> = -ra*b3>
p_dao_dbo> =-l*da3>
p_dbo_cb> = rb*unitvec(a3> - dot(a3>, db2>)*db2>)

f = dot(p_ca_cb>, a3>)
dfdq1 = d(f, q1)
dfdq2 = d(f, q2)

output f, dfdq1, dfdq2
code algebraic() slotted_discs_hc_al.c
