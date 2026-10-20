"""Support-4 system, complex side, the 16 coordinate totals of the Gram decoder: for ANY F_2 labels the cost of a
total is at least the size of a pairwise non-orthogonal set inside its read set (those labels are independent).
Also the cost with the labels N = I + A_orth (h = 111)."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linesys
from survey import clique
L = linesys.rm24w4(); v = len(L); G = L @ L.T
def f2rank(M):
    P = np.packbits((M & 1).astype(np.uint8), axis=1); basis = {}
    for r in P:
        x = int.from_bytes(r.tobytes(), 'big')
        while x:
            hb = x.bit_length() - 1
            if hb in basis: x ^= basis[hb]
            else: basis[hb] = x; break
    return len(basis)
N = np.eye(v, dtype=np.int64) + (G == 0).astype(np.int64)
lb, cost = [], []
for a in range(16):
    R = np.nonzero(L[:, a])[0]
    lb.append(len(clique(L[R], seeds=(1, 2), restarts=15)))
    cost.append(f2rank(N[R, :]))
res = dict(read_set_size=int(np.count_nonzero(L[:, 0])), clique_lower_bound_per_total=lb, sum_lower=sum(lb), x_lower=sum(lb) / v,
           cost_with_I_plus_Aorth=cost, x_with_I_plus_Aorth=sum(cost) / v)
print(res); json.dump(res, open('rm24w4_coord.json', 'w'), indent=1)
