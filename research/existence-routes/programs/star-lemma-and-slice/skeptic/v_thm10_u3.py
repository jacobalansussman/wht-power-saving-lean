# v_thm10_u3.py (v-e3-slice6): Theorem 10 on the genuine (numerical) U_3 unit, and the local chain at block 3.
import numpy as np, itertools
d = np.load('../u3/u3_solution.npz'); A, B, G = d['A'], d['B'], d['Gam']; dirs = [int(z) for z in d['dirs']]; W = 2; m = 3
E = np.eye(7) + G - B @ A
mons = [(j, k) for j in range(m) for k in range(j, m)]
worst_even = 1e9; worst_odd = 0
for q in itertools.product([0, 1], repeat=6):
    D = [i for i, z in enumerate(dirs) if sum(q[t] * ((z >> j) & 1) * ((z >> k) & 1) for t, (j, k) in enumerate(mons)) % 2]
    tr = sum(q[t] for t, (j, k) in enumerate(mons) if j == k) % 2
    if not D: continue
    sv = np.linalg.svd(E[np.ix_(D, D)], compute_uv=False)
    if tr == 0: worst_even = min(worst_even, sv.min())
    else: worst_odd = max(worst_odd, sv[len(D) - W:].max() if len(D) > W else 0)
print("Theorem 10 on U_3: min smallest singular value over even forms %.3g (nonsingular); max singular value beyond slack on odd forms %.2g" % (worst_even, worst_odd))
Fi = [i for i, z in enumerate(dirs) if z & 1 and z & 2]
print("F = {z1=z2=1} moves", [dirs[i] for i in Fi], "; E[F,F] =", np.round(E[np.ix_(Fi, Fi)], 6).tolist())
