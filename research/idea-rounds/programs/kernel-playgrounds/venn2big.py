# float-rank screen of Lambda^a (x) Lambda^b Venn certificates for the larger (a,b) at block 3; best ones re-ranked mod p.
import sys, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hookk import rank_mod, wedge, kron_list
P = 2147483647
def frank(B):
    if B.size == 0: return 0
    return int(np.linalg.matrix_rank(B.astype(np.float64)))
m = 3; n = 8
KS = {'H': [[1,1],[1,-1]], 'L': [[1,0],[1,1]]}
pairs = [(int(a), int(b)) for a, b in (s.split(',') for s in sys.argv[1:])]
for name, K in KS.items():
    K = np.array(K, dtype=object); I = np.eye(2, dtype=object)
    K1 = kron_list([K, I, I]); T = kron_list([K, K, K])
    for (a, b) in pairs:
        SA = list(itertools.combinations(range(n), a)); SB = list(itertools.combinations(range(n), b))
        inter = np.array([len(set(A) & set(B)) for A in SA for B in SB])
        vals = sorted(set(inter.tolist()))
        Ms = {}
        for nm, g in (('cap', K1), ('dem', T)):
            gi = np.array([[int(v) for v in row] for row in g], dtype=np.int64)
            gp = gi % P
            # exact integer wedge via mod-P minors mapped back to symmetric range (entries are small)
            wa = wedge(gp, a, P, SA); wb = wedge(gp, b, P, SB)
            wa = np.where(wa > P // 2, wa - P, wa); wb = np.where(wb > P // 2, wb - P, wb)
            Ms[nm] = np.kron(wa, wb)
        best = None
        for r in range(1, len(vals)):
            for S in itertools.combinations(vals, r):
                inY = np.isin(inter, S)
                cap = frank(Ms['cap'][np.ix_(~inY, inY)]); dem = frank(Ms['dem'][np.ix_(~inY, inY)])
                rat = dem / cap if cap else 0
                if best is None or rat > best[0]: best = (rat, S, cap, dem, int(inY.sum()))
        print(f"{name} a={a} b={b} dim={len(inter)}: best float ratio {best[0]:.4f} S={best[1]} cap={best[2]} dem={best[3]} dimY={best[4]}", flush=True)
