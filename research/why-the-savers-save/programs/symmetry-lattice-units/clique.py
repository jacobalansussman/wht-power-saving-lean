"""Lower bound for the F_2 minrank of the BW16 orthogonality pattern: a large set of pairwise NON-orthogonal lines
(its principal submatrix is forced to be the identity).  Randomised greedy + simple local search (plateau moves)."""
import numpy as np, random, sys
L = np.load('bw16_L.npy'); v = len(L)
G = (L @ L.T) // 2
NO = (G != 0); np.fill_diagonal(NO, False)
nb = [set(np.nonzero(NO[i])[0].tolist()) for i in range(v)]
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
best = []
def greedy(start):
    C = [start]; cand = set(nb[start])
    while cand:
        # pick candidate with most neighbours inside cand (classic heuristic), random tie-break
        sc = max(cand, key=lambda u: (len(nb[u] & cand), rng.random()))
        C.append(sc); cand &= nb[sc]
    return C
for it in range(60):
    C = greedy(rng.randrange(v))
    if len(C) > len(best): best = C; print('clique', len(best), flush=True)
# local search: (1,1)-swaps then extend
Cs = set(best)
for it in range(3000):
    # vertices adjacent to all but one of Cs
    out = [u for u in range(v) if u not in Cs and len(nb[u] & Cs) == len(Cs) - 1]
    if not out: break
    u = rng.choice(out); w = next(x for x in Cs if x not in nb[u])
    Cs.remove(w); Cs.add(u)
    add = [x for x in range(v) if x not in Cs and Cs <= nb[x]]
    while add:
        x = rng.choice(add); Cs.add(x); add = [y for y in add if y != x and y in nb[x]]
    if len(Cs) > len(best): best = list(Cs); print('clique (ls)', len(best), flush=True)
assert all(G[a, b] != 0 for a in best for b in best)
print('best pairwise non-orthogonal set:', len(best))
np.save('clique_best.npy', np.array(sorted(best)))
