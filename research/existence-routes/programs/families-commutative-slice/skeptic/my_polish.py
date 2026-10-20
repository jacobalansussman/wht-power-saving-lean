import numpy as np
exec(open("my_solve_u3.py").read().split("rng = np.random")[0])
d = np.load("u3_solution.npz")
A, B, G = d["A"], d["B"], d["Gam"]
v = np.concatenate([A.ravel(), B.ravel(), np.array([G[j, i] for (j, i) in low])])
v = np.concatenate([v.real, v.imag])
print("start residual %.3e" % np.linalg.norm(res(v)))
for it in range(3):
    s = least_squares(res, v, method="trf", xtol=3e-16, ftol=3e-16, gtol=3e-16, max_nfev=3000)
    v = s.x
    print("residual %.3e" % np.linalg.norm(s.fun))
A, B, G = unpack(v)
np.savez("u3_solution.npz", A=A, B=B, Gam=G, dirs=np.array(dirs), W=W, m=m)
