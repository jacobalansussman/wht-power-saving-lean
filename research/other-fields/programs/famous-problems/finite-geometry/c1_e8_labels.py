# c1: the community's "E8" label family (README r6 section 6) as a piece of the E8-mod-2 quadric.
# Claim [derived here]: labels = odd-weight words of F_2^9 of weight 3 or 7; "know" = even overlap.
# Map x -> x + 1 (complement) sends labels onto ALL vectors of the even-weight space E (dim 8) with
# q(x) = wt(x)/2 mod 2 = 1 (the 120 nonsingular points of O+(8,2) = E8/2E8 roots mod sign), and
# "know" (x.y = 0) becomes "non-orthogonal" (B = 1). Checks: counts, quadric type (135 singular),
# the flip, SRG parameters, and the same graph built from the 240 E8 roots.
from itertools import combinations, product
n = 9
def wt(x): return bin(x).count("1")
def dot(x, y): return wt(x & y) & 1
ONE = (1 << n) - 1
labels = [x for x in range(1 << n) if wt(x) in (3, 7)]
print("labels", len(labels), "(weight 3:", sum(wt(x) == 3 for x in labels), ", weight 7:", sum(wt(x) == 7 for x in labels), ")")
know = {x: {y for y in labels if y != x and dot(x, y) == 0} for x in labels}
print("knowing degrees", sorted({len(v) for v in know.values()}))
E = [x for x in range(1 << n) if wt(x) % 2 == 0]
q = lambda x: (wt(x) // 2) & 1
print("dim E = 8:", len(E) == 256, "; singular nonzero in E:", sum(1 for x in E if x and q(x) == 0), "(O+(8,2) has 135, O-(8,2) has 119)")
nons = sorted(x for x in E if q(x) == 1)
img = sorted(x ^ ONE for x in labels)
print("labels^1 == all nonsingular points of E:", img == nons, len(nons))
print("know <=> images non-orthogonal:", all((dot(x, y) == 0) == (dot(x ^ ONE, y ^ ONE) == 1) for x, y in combinations(labels, 2)))
def srg(V, adj):
    k = {len(adj[v]) for v in V}
    lam = {len(adj[a] & adj[b]) for a, b in combinations(V, 2) if b in adj[a]}
    mu = {len(adj[a] & adj[b]) for a, b in combinations(V, 2) if b not in adj[a]}
    return len(V), k, lam, mu
print("knowing graph (v,k,lambda,mu):", srg(labels, know))
# E8 roots: 112 of type (+-1,+-1,0^6), 128 of type (+-1/2)^8 with an even number of minus signs. Work in doubled coords.
roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product((2, -2), repeat=2):
        r = [0] * 8; r[i] = si; r[j] = sj; roots.append(tuple(r))
for s in product((1, -1), repeat=8):
    if sum(1 for t in s if t < 0) % 2 == 0: roots.append(s)
print("E8 roots", len(roots))
lines = sorted({max(r, tuple(-t for t in r)) for r in roots})
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4   # true inner product
adjE = {a: {b for b in lines if b != a and ip(a, b) != 0} for a in lines}
print("E8 line graph, adjacency = non-orthogonal (v,k,lambda,mu):", srg(lines, adjE))
