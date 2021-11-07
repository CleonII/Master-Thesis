function ODE_model!(du, u, p, t)
k1 = p
c1 = 2.0
s1, s1s2, s2 = u
du[1] = -1.0 * (c1 * k1 * s1 * s2)
du[2] = +1.0 * (c1 * k1 * s1 * s2)
du[3] = -1.0 * (c1 * k1 * s1 * s2)
nothing
end
