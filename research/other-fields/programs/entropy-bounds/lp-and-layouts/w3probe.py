# w3probe.py -- look at the cuts the optimal packings use (block 2/3): are they level sets {dim(L cap X) >= j}?
import sys, random, numpy as np, itertools
import w3cut, dpg
from fractions import Fraction
m = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3])
g = w3cut.G(m); rng = random.Random(seed)
N = 1 << (2 * m)
# all subspaces of F_2^(2m) as frozensets
subs = set()
def span(vs):
    s = {0}
    for v in vs: s |= {x ^ v for x in s}
    return frozenset(s)
frontier = {frozenset({0})}; subs |= frontier
while frontier:
    new = set()
    for S in frontier:
        for v in range(1, N):
            if v not in S:
                T = frozenset(S | {x ^ v for x in S})
                if T not in subs: new.add(T)
    subs |= new; frontier = new
subs = [S for S in subs if 1 < len(S) < N]
print(len(subs), "proper nonzero subspaces")
dims = {}
for S in subs:
    dims[S] = np.array([len(S & L).bit_length() - 1 for L in g.verts])
def name(inQ):
    for side, Q in (('', inQ), ('~', ~inQ)):
        for S in subs:
            d = dims[S]
            for j in range(1, len(S).bit_length()):
                if np.array_equal(d >= j, Q):
                    iso = all(dpg.form(a, b, m) == 0 for a in S for b in S)
                    return f"{side}[dim(L^X)>={j}, dimX={len(S).bit_length()-1}{' iso' if iso else ''}]"
    return f"other(|Q|={int(inQ.sum())})"
for tr in range(trials):
    pairs = w3cut.rand_3k2(g, rng)
    val, lam, Qs, rhs, rnd = w3cut.full_cut_lp(g, pairs)
    sup = [i for i in range(len(lam)) if lam[i] > 1e-9]
    print(f"trial {tr}: {pairs} value {val:.6f}")
    for i in sup:
        print(f"   lam {Fraction(lam[i]).limit_denominator(1000)}  sep {rhs[i]}  {name(Qs[i])}")
