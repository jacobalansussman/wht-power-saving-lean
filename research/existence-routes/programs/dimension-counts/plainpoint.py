# Jacobian rank of the program map AT a plain schedule (a known point of the unit variety), against the generic rank.
# Block 2, rep VV (and VL2): plain 1-in-2:  u1 = D E,  u2 = D (pi^-1 u1),  O = pi u2 = K E   (pi swaps the two bits).
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dimcount import Shape, rank_sv, letters, perm_mat
rep = sys.argv[1] if len(sys.argv) > 1 else 'VV'
S = Shape(2, rep, 'engine', 1, [1, 1], [2], None, 0)
f = __import__('dimcount').rep_fn(rep)
pi = f(perm_mat(4, [(x >> 1) | ((x & 1) << 1) for x in range(4)]))
D = S.Wr[1]; K = S.Wr[2]
assert np.allclose(pi @ D @ np.linalg.inv(pi) @ D, K) or np.allclose(np.linalg.inv(pi) @ D @ pi @ D, K)
def coords(Mtx):            # coordinates of an element of A in the (orthonormal Frobenius) basis
    return np.einsum('aij,ij->a', S.A.conj(), Mtx)
Z = np.zeros(S.s, complex)
# slots: step1 [E], step2 [E, u1], out [E, u1, u2]
coef = [np.array([coords(np.eye(S.d))]), np.array([Z, coords(np.linalg.inv(pi))]), np.array([Z, Z, coords(pi)])]
gens, outs = S.forward(coef)
print('residual at plain point:', float(np.abs(outs[0] - S.target(0)).max()))
J = S.jac_sketch(coef, S.W * S.W * S.gU + 24)
r0, gap0 = rank_sv(J)
Jg = S.jac_sketch(S.random_point(), S.W * S.W * S.gU + 24)
rg, gapg = rank_sv(Jg)
print(json.dumps({'rep': rep, 'shape': 'plain 1-in-2', 'P': S.P, 'Eeff': S.W * S.W * S.gU, 'rank_at_plain_point': r0, 'gap': gap0,
                  'generic_rank': rg, 'gap_generic': gapg, 'local_dim_if_smooth': S.P - r0}))
