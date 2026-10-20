# Independent enumeration of the block-3 frame graph: Lagrangians of F_2^6 with omega((x|z),(x'|z')) = x.z' + x'.z.
# Counts frames, edges (dim intersect = 2), intervals (opposite pairs), and per-edge end / middle incidences.
from itertools import combinations
def om(a, b):
    x, z = a >> 3, a & 7; x2, z2 = b >> 3, b & 7
    return (bin(x & z2).count("1") + bin(x2 & z).count("1")) & 1
def span(vs):
    S = {0}
    for v in vs:
        if v not in S: S |= {s ^ v for s in S}
    return frozenset(S)
L = set()
for a, b, c in combinations(range(1, 64), 3):
    if om(a, b) or om(a, c) or om(b, c): continue
    S = span([a, b, c])
    if len(S) == 8: L.add(S)
L = sorted(L, key=lambda s: sorted(s)); n = len(L); idx = {s: i for i, s in enumerate(L)}
dim = lambda S: len(S).bit_length() - 1
def d(i, j): return 3 - dim(L[i] & L[j])
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if d(i, j) == 1]
eid = {e: k for k, e in enumerate(edges)}
opp = [(i, j) for i in range(n) for j in range(i + 1, n) if d(i, j) == 3]
print("frames", n, "edges", len(edges), "opposite pairs (intervals)", len(opp))
end = [0] * len(edges); mid = [0] * len(edges); sizes = set()
for (a, b) in opp:
    I = [k for k in range(n) if d(a, k) + d(k, b) == 3]
    lev = {k: d(a, k) for k in I}
    sizes.add(len(I))
    E_I = [(min(i, j), max(i, j)) for i in I for j in I if i < j and d(i, j) == 1]
    for e in E_I:
        l = sorted((lev[e[0]], lev[e[1]]))
        if l == [0, 1] or l == [2, 3]: end[eid[e]] += 1
        elif l == [1, 2]: mid[eid[e]] += 1
        else: raise SystemExit("non-covering edge in interval")
    assert len(E_I) == 35, len(E_I)
print("interval sizes", sizes, "end incidences per edge", set(end), "middle incidences per edge", set(mid))
# G3 arithmetic: sum over intervals of (mu + k/2) >= 2 per interval -> (mid + end/2) E >= 2 * #intervals
m0, e0 = mid[0], end[0]
print("bound E >= 2*%d/(%d + %d/2) = %.4f" % (len(opp), m0, e0, 2 * len(opp) / (m0 + e0 / 2)))
