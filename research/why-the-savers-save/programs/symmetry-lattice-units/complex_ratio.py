"""Complex-type BW16 unit (F_2 labels, rational decoder B = Gram, the fat-label 'ratio' = cst/v of the E8 0.600 /
E7 0.778 table).  Labels: N = I + A (A = all orthogonal pairs, h = rank 191) or N = frame blocks (h = 135).
Decoder: the 16 coordinate totals (Gram = sum over coordinates); cost of a total = F_2 rank of the labels it reads
= rank N[R, :].  Also the pairwise non-orthogonal set (clique) that forces h >= 72 for every label choice."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import linesys
L = linesys.bw16(); v = len(L); G = L @ L.T
def f2rank(M):
    P = np.packbits((M & 1).astype(np.uint8), axis=1)
    basis = {}
    for r in P:
        x = int.from_bytes(r.tobytes(), 'big')
        while x:
            hb = x.bit_length() - 1
            if hb in basis: x ^= basis[hb]
            else: basis[hb] = x; break
    return len(basis)
O = (G == 0).astype(np.int64); I = np.eye(v, dtype=np.int64)
# the 15 'special' orthogonal partners of each line (common |ip|=2 count 56) close up into 135 frames of 16
A2 = (np.abs(G) == 4).astype(np.int64)             # |ip| = 2 at the norm-4 scale (norm-8 ips are 0, +-2, +-4)
vals = sorted(set(np.abs(G[np.triu_indices(v, 1)]).tolist()))
res = dict(abs_ips_norm8=vals)
labels = {'I+A_orth': I + O}
# frames: connected components of the special-orthogonal relation = lines whose difference/sum is divisible by 2 (x = y mod 2)
C2 = A2 @ A2
spec = (G == 0) & (C2 == 56); np.fill_diagonal(spec, False)
deg = spec.sum(1); assert set(deg.tolist()) == {15}
Nf = spec.astype(np.int64) + I
comps = {tuple(np.nonzero(Nf[i])[0].tolist()) for i in range(v)}
assert len(comps) == 135 and all(len(c) == 16 for c in comps)
res['frames_count'] = len(comps)
labels['frames'] = Nf
for name, N in labels.items():
    assert np.all(np.diag(N) == 1) and not np.any(N[(G != 0) & ~np.eye(v, dtype=bool)] % 2), name
    h = f2rank(N)
    costs = [f2rank(N[np.nonzero(L[:, a])[0], :]) for a in range(16)]
    res[name] = dict(h=h, totals=16, cost_per_total=costs, cst=sum(costs), x=sum(costs) / v)
    print(name, res[name], flush=True)
json.dump(res, open('complex_ratio.json', 'w'), indent=1)
