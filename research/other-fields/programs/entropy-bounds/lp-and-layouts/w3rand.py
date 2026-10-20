# w3rand.py -- do three-pair min-cost cut gaps exist in general graphs? (random small graphs, random 3K2 demands)
import sys, random, numpy as np, itertools, time
from scipy.sparse.csgraph import shortest_path
import scipy.sparse as sp
import w3image
n = int(sys.argv[1]); p = float(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4])
rng = random.Random(seed); t0 = time.time(); gaps = 0
for tr in range(trials):
    E = [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < p]
    A = sp.csr_matrix((np.ones(len(E)), ([u for u, v in E], [v for u, v in E])), shape=(n, n))
    D = shortest_path(A, directed=False, unweighted=True)
    if np.isinf(D).any(): continue
    vs = rng.sample(range(n), 6); P = [(vs[0], vs[1]), (vs[2], vs[3]), (vs[4], vs[5])]
    route = sum(D[s, t] for s, t in P)
    val = w3image.cut_value(n, E, P)[0]
    if val < route - 1e-6:
        gaps += 1; print(f"GAP trial {tr}: n={n} |E|={len(E)} P={P} route {route} cut {val:.4f}  E={E}", flush=True)
print(f"done: {gaps} gaps in {trials} trials ({time.time()-t0:.0f}s)")
