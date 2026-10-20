# w3rand2.py -- three-pair min-cost cut gaps in general graphs when every pair is DIAMETRAL (d(s_i,t_i) = diameter),
# as opposite frames are in Gamma_m. Also report whether some two pairs are compatible.
import sys, random, numpy as np, itertools, time
from scipy.sparse.csgraph import shortest_path
import scipy.sparse as sp
import w3image
n = int(sys.argv[1]); p = float(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed); t0 = time.time(); gaps = 0; tested = 0; incompat = 0
def compat(D, m, P):
    for i, j in ((0, 1), (0, 2), (1, 2)):
        (si, ti), (sj, tj) = P[i], P[j]
        if (D[si, tj] == m and D[sj, ti] == m) or (D[si, sj] == m and D[ti, tj] == m): return True
    return False
for tr in range(trials):
    E = [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < p]
    A = sp.csr_matrix((np.ones(len(E)), ([u for u, v in E], [v for u, v in E])), shape=(n, n))
    D = shortest_path(A, directed=False, unweighted=True)
    if np.isinf(D).any(): continue
    m = D.max()
    diam = [(u, v) for u in range(n) for v in range(u + 1, n) if D[u, v] == m]
    # three vertex-disjoint diametral pairs
    for _ in range(50):
        P = rng.sample(diam, 3) if len(diam) >= 3 else None
        if P and len({x for q in P for x in q}) == 6: break
    else: continue
    tested += 1
    c = compat(D, m, P)
    if not c: incompat += 1
    val = w3image.cut_value(n, E, P)[0]
    if val < 3 * m - 1e-6:
        gaps += 1; print(f"GAP: n={n} diam {m} P={P} cut {val:.4f} compatible={c} E={E}", flush=True)
print(f"done: {tested} tested, {incompat} with no compatible pair of pairs, {gaps} gaps ({time.time()-t0:.0f}s)")
