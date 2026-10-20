"""base.py (d-isocut) -- P112 ports, F_2^9 labels, orthogonality, octads, totally isotropic (t.i.) cuts.
Ports are lines 8..119 of g-break/common.py (lines 0..7 = removed octad).  S,T interfere (orthogonal roots) iff
u_S.u_T = 1 (S != T); S knows T iff u_S.u_T = 0."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "../copy-wall-idea-probes/break"))
from common import *          # L, G, B0, U, UI, f2span, f2dim, dot9
P = list(range(8, 120))
IDX = {S: i for i, S in enumerate(P)}
lab = [UI[S] for S in P]                       # labels of the 112 ports, index 0..111
N = len(P)
orth = [set(k for k in range(N) if k != i and G[P[i]][P[k]] == 0) for i in range(N)]
for i in range(N):
    for k in orth[i]:
        assert dot9(lab[i], lab[k]) == 1
    for k in range(N):
        if k != i and k not in orth[i]:
            assert dot9(lab[i], lab[k]) == 0
    assert dot9(lab[i], lab[i]) == 1

def octads():
    out = set()
    def rec(cur, cand):
        if len(cur) == 8:
            out.add(tuple(cur)); return
        for T in sorted(cand):
            if T < cur[-1]: continue
            rec(cur + [T], cand & orth[T])
    for S in range(N):
        rec([S], set(T for T in orth[S] if T > S))
    return sorted(out)

def omega(Vs):
    Vs = set(Vs); best = [0]
    def r(size, cand):
        if size + len(cand) <= best[0]: return
        if not cand:
            best[0] = max(best[0], size); return
        T = min(cand)
        r(size + 1, cand & orth[T]); r(size, cand - {T})
    r(0, Vs); return best[0]

def span_ints(basis):
    s = {0}
    for b in basis:
        s |= {x ^ b for x in s}
    return s

def ti_subspaces(maxdim=4):
    """all totally isotropic subspaces (as frozensets of nonzero vectors + a canonical basis), by dim"""
    iso = [h for h in range(1, 512) if not dot9(h, h)]
    res = {1: {frozenset([h]): [h] for h in iso}}
    for d in range(2, maxdim + 1):
        res[d] = {}
        for key, basis in res[d - 1].items():
            sp = span_ints(basis)
            for h in iso:
                if h in sp: continue
                if any(dot9(h, b) for b in basis): continue
                nb = basis + [h]
                k2 = frozenset(span_ints(nb) - {0})
                if k2 not in res[d]:
                    res[d][k2] = nb
    return res

def C_of(basis):
    return [i for i in range(N) if all(not dot9(lab[i], b) for b in basis)]
