autoz off

newtonian n
bodies s
frames da, db, dagl, dbgl
frames a, b, cam

constants m, ra, rb, l, g, alpha
constants k

variables q{5}'
variables w'
variables fax, fay, faz, fbx, fby, fbz

points ca, cb

simprot(n, a, 3, q1)
simprot(a, b, 1, q2)
simprot(b, da, 2, q3)
simprot(da, db, 3, alpha)

p_ca_dao> = express(-ra*b3>, a)
p_dao_dbo> = express(-l*da3>, a)
p_dbo_cb> = express(rb*unitvec(a3> - dot(a3>, db2>)*db2>), a)

f = dot(p_ca_cb>, a3>)
dfdq2 = d(f, q2)
dfdq3 = d(f, q3)

%constants r
%hc = evaluate(dot(p_ca_cb>, a3>), alpha=pi/2, l=r, ra=r, rb=r)

output f, dfdq2, dfdq3
code algebraic() slotted_discs_hc_al.c
