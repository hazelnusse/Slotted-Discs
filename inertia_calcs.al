bodies a, b
frames s
constants ma, ra, mb, rb, l, alpha, Ia, Ja, Ib, Jb

simprot(a, b, 1, alpha)

inertia a, Ia, Ia, Ja, 0, 0, 0
inertia b, Ib, Ib, Jb, 0, 0, 0
mass a=ma
mass b=mb

p_ao_bo> = -l*a1>
p_ao_so> = cm(ao)

k = dot(p_ao_so>, a1>)

I_S_SO>> = I_A_SO>> + I_B_SO>>

Ixx = dot(a1>, dot(I_S_SO>>, a1>))
Iyy = dot(a2>, dot(I_S_SO>>, a2>))
Izz = dot(a3>, dot(I_S_SO>>, a3>))
Iyz = dot(a2>, dot(I_S_SO>>, a3>))
Izy = dot(a3>, dot(I_S_SO>>, a2>))
output Ixx, Iyy, Izz, Iyz, k
code algebraic() inertias_al.c
