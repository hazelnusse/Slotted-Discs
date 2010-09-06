bodies da, db
frames s
constants ma, ra, mb, rb, l, alpha, Ia, Ja, Ib, Jb

simprot(da, db, 1, alpha)
simprot(da, s, 1, 0)

inertia da, Ia, Ia, Ja, 0, 0, 0
inertia db, Ib, Ib, Jb, 0, 0, 0
mass da=ma
mass db=mb

p_dao_dbo> = -l*s1>
p_dao_so> = cm(dao, da, db)
k = -dot(p_dao_so>, s1>)

I_S_SO>> = I_DA_SO>> + I_DB_SO>>

Ixx = dot(s1>, dot(I_S_SO>>, s1>))
Iyy = dot(s2>, dot(I_S_SO>>, s2>))
Izz = dot(s3>, dot(I_S_SO>>, s3>))
Iyz = dot(s2>, dot(I_S_SO>>, s3>))
Izy = dot(s3>, dot(I_S_SO>>, s2>))
output Ixx, Iyy, Izz, Iyz, k
code algebraic() inertias_al.c
