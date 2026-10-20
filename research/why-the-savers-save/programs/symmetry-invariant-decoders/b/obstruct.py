"""obstruct.py (t-km): exact combinatorial obstructions for a frame family on P112 (any payload field).
J(T) = frames containing port T; rho_T in F^K is supported on J(T) and nonzero; fit needs rho_S not in
span{rho_T : T orthogonal to S}.  Tests:
 (i)  two orthogonal ports with J(S) = J(T) = {j}  -> rho_S, rho_T parallel -> impossible;
 (ii) clique Hall on every maximal clique of the orthogonality graph inside P (not only octads): rho on a clique
      must be independent, so its ports need a matching into distinct frames;
 (iii) for each S: if J(S) is inside the union of supports of single-frame-forced orthogonal ports... (span test
      with 'forced' vectors: ports with |J| = 1 have rho = e_j up to scale; S fails if e_{J(S)} all forced).
usage: python -I obstruct.py fam63.json"""
import sys, json, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from frames import *
import numpy as np
fams = json.load(open(sys.argv[1]))
OM = [sum(1 << T for T in range(NP) if ORTH[S, T]) for S in range(NP)]
def maxcliques():
    out = []
    def bk(R, Pm, Xm):
        if not Pm and not Xm: out.append(R); return
        u = (Pm | Xm).bit_length() - 1
        c = Pm & ~OM[u]
        while c:
            low = c & -c; v = low.bit_length() - 1; c ^= low
            bk(R + [v], Pm & OM[v], Xm & OM[v]); Pm &= ~low; Xm |= low
    bk([], (1 << NP) - 1, 0)
    return out
MC = maxcliques()
from collections import Counter
print('maximal cliques of the orthogonality graph on P112:', len(MC), Counter(len(c) for c in MC))
for f in fams:
    S = np.zeros((len(f['frames']), NP), bool)
    for k, fr in enumerate(f['frames']): S[k, fr] = True
    J = [frozenset(np.nonzero(S[:, T])[0].tolist()) for T in range(NP)]
    single = [T for T in range(NP) if len(J[T]) == 1]
    bad1 = [(a, b) for a, b in itertools.combinations(single, 2) if ORTH[a, b] and J[a] == J[b]]
    print('family dims', f['dims'], ': ports with one frame', len(single), '; orthogonal pairs sharing their only frame:', len(bad1), bad1[:5])
    # clique Hall on all maximal cliques
    def hall_ok(cl):
        # bipartite matching cl -> frames
        match = {}
        def aug(T, seen):
            for j in J[T]:
                if j in seen: continue
                seen.add(j)
                if j not in match or aug(match[j], seen): match[j] = T; return True
            return False
        return all(aug(T, set()) for T in cl)
    fails = [c for c in MC if not hall_ok(c)]
    print('   maximal cliques failing Hall:', len(fails), 'of', len(MC), [len(c) for c in fails[:10]])
    # forced test (iii): rho_T = e_j for single-frame T; S fails if every j in J(S) is hit by an orthogonal single-frame T
    bad3 = [Sx for Sx in range(NP) if all(any(J[T] == {j} for T in range(NP) if ORTH[Sx, T]) for j in J[Sx])]
    print('   ports whose every frame is pinned by an orthogonal single-frame port:', len(bad3), bad3[:10])
