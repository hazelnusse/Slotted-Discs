bodies da, db
frames s
constants ma, ra, mb, rb, l, alpha, Ia, Ja, Ib, Jb

simprot(da, db, 3, alpha)
da_s = [1, 0, 0; 0, 1, 0; 0, 0, 1]

inertia da, Ia, Ja, Ia, 0, 0, 0
inertia db, Ib, Jb, Ib, 0, 0, 0
mass da=ma
mass db=mb

p_dao_dbo> = -l*da3>
p_dao_so> = cm(dao, da, db)
k = -dot(p_dao_so>, da3>)

I_S_SO>> = I_DA_SO>> + I_DB_SO>>

Ixx = dot(s1>, dot(I_S_SO>>, s1>))
Iyy = dot(s2>, dot(I_S_SO>>, s2>))
Izz = dot(s3>, dot(I_S_SO>>, s3>))
Ixy = dot(s1>, dot(I_S_SO>>, s2>))
output Ixx, Iyy, Izz, Ixy, k
code algebraic() inertias_al.c
