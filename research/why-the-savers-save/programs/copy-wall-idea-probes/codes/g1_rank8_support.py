"""g1 (g-codes): can ANY field (Q, F_3, F_4, ..., or vector-linear via F_{2^t}) give a rank-8 decoder on the
BIT E8 pattern (B[S,T]=0 when roots non-orthogonal, B[S,S]!=0)?  Field-free support propagation.
Gauge: K = 8 pairwise non-orthogonal lines -> B[K,K] diagonal, X_K = I, Y_K = D.  Then X_S, Y_S are supported on
Z(S) = {k in K : S orth k}.  Forced zeros propagate:  |Z(S)cap Z(T)| = 1 = {k} for a non-orth pair S,T gives
X_S[k] Y_T[k] = 0.  If a line ends with no admissible k with X_S[k],Y_S[k] both possibly nonzero -> UNSAT for
every field.  Run: python3 -I g1_rank8_support.py"""
import sys, itertools
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
from copyfloor import e8_lines
L, G = e8_lines()
n = 120
nonorth = [[S != T and G[S][T] != 0 for T in range(n)] for S in range(n)]
# all maximal cliques of pairwise non-orthogonal lines, size 8 (and check 9)
adj = [sum(1 << T for T in range(n) if nonorth[S][T]) for S in range(n)]
def cliques(size):
    out = []
    def rec(cur, cand):
        if len(cur) == size:
            out.append(tuple(cur)); return
        c = cand
        while c:
            low = c & -c; i = low.bit_length() - 1; c ^= low
            rec(cur + [i], cand & adj[i] & ~((low << 1) - 1))
            if len(out) > 200000: return
    rec([], (1 << n) - 1)
    return out
c9 = cliques(9)
print('cliques of 9 pairwise non-orthogonal lines:', len(c9))
c8 = cliques(8)
print('cliques of 8 (capped 200000):', len(c8))
def test(K):
    Ks = set(K)
    Z = {S: [k for k in K if G[S][k] == 0] for S in range(n) if S not in Ks}
    # possible[S] = set of k with X_S[k] possibly nonzero; possY same
    pX = {S: set(Z[S]) for S in Z}; pY = {S: set(Z[S]) for S in Z}
    for k in K:
        pX[k] = {k}; pY[k] = {k}
    allS = list(range(n))
    changed = True
    while changed:
        changed = False
        for S in allS:
            # forced nonzero: if exactly one k admissible for both X and Y, X_S[k] and Y_S[k] are nonzero
            both = pX[S] & pY[S]
            if not both:
                return ('UNSAT', S)
            if len(both) == 1 and len(pX[S]) >= 1:
                (k,) = both
                # X_S[k] != 0 and Y_S[k] != 0 (the diagonal sum has only that term)
                for T in allS:
                    if nonorth[S][T]:
                        common = (pX[S] & pY[T])
                        # constraint sum_{k' in common} X_S[k'] Y_T[k'] = 0; if common == {k}: Y_T[k] = 0
                        if common == {k} and k in pY[T]:
                            pY[T].discard(k); changed = True
                        common2 = (pX[T] & pY[S])
                        if common2 == {k} and k in pX[T]:
                            pX[T].discard(k); changed = True
    sizes = sorted(len(pX[S] & pY[S]) for S in allS)
    return ('open', sizes[:10], sum(sizes))
from collections import Counter
res = Counter()
for i, K in enumerate(c8[:2000]):
    r = test(K)
    res[r[0]] += 1
    if i < 3: print(K, r)
print('bit E8 pattern, rank 8, support propagation over 2000 gauge cliques:', dict(res))
# Z-size profile for the first clique
K = c8[0]
print('Z-size histogram (lines outside K):', Counter(len([k for k in K if G[S][k] == 0]) for S in range(n) if S not in K))
