"""The support-4 Clifford system (140 weight-4 words of RM(2,4) x 8 sign classes = 1120 lines in Q^16): the one system
whose clique bound (h >= 40) does not by itself kill the complex side at s = 3.  Upper bounds for its F_2 minrank from
invariant matrices, a longer clique search, and the copy charge a complex unit with that window would need."""
import sys, os, json, random, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linesys
from fractions import Fraction as Fr
from ceil_grid import Xf, x_needed
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
I = np.eye(v, dtype=np.int64); O = (G == 0).astype(np.int64)
nonorth = (G != 0) & ~np.eye(v, dtype=bool)
res = dict(v=v, ips=sorted(set(np.abs(G[np.triu_indices(v, 1)]).tolist())))
res['rank_I_plus_Aorth'] = f2rank(I + O)
# refine orthogonal pairs by common-neighbour counts and try every union of classes
A = nonorth.astype(np.int64); C = A @ A
cls = sorted(set(C[O == 1].tolist())); res['orth_classes_by_common_nonorth'] = cls
best = None
for mask in range(1, 1 << len(cls)):
    keep = np.isin(C, [c for i, c in enumerate(cls) if mask >> i & 1]) & (O == 1)
    N = I + keep.astype(np.int64); r = f2rank(N)
    if best is None or r < best[0]: best = (r, [c for i, c in enumerate(cls) if mask >> i & 1])
res['best_invariant_minrank_upper'] = best
cq = clique(L, seeds=tuple(range(1, 11)), restarts=40); res['clique'] = len(cq)
h_lo, h_hi = len(cq), best[0]
for h in (h_lo, h_hi):
    res['x_needed_h%d' % h] = {s: x_needed(h, Xf(s, 9, Fr(64, 112)), s=s) for s in (3, 5)}
    res['x_needed_h%d_vs_E8best' % h] = x_needed(h, 5.186252252886545e-3)
print(json.dumps(res, indent=1))
json.dump(res, open('rm24w4.json', 'w'), indent=1)
