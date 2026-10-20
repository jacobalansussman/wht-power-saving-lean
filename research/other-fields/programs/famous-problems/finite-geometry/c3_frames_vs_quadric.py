# c3: the 135 block-3 frames (Lagrangian subspaces of symplectic F_2^6 = dual polar space DSp(6,2))
# against the 135 singular points of the E8-mod-2 quadric O+(8,2) (same space as c1).
# Compare the "opposite frames" graph with the "non-orthogonal singular points" graph (SRG parameters),
# and the "frame distance 1" graph with nothing (valency 14 has no analogue on the quadric).
from itertools import combinations
def span(vs):
    S = {0}
    for v in vs: S |= {s ^ v for s in S}
    return frozenset(S)
def symp(a, b):   # F_2^6 = (x1..x3 | z1..z3), omega = x.z' + z.x'
    return (bin((a & 7) & (b >> 3)).count("1") + bin((a >> 3) & (b & 7)).count("1")) & 1
L = set()
vecs = range(1, 64)
for a, b, c in combinations(vecs, 3):
    if symp(a, b) or symp(a, c) or symp(b, c): continue
    S = span([a, b, c])
    if len(S) == 8: L.add(S)
L = sorted(L, key=sorted)
print("Lagrangians (frames at block 3):", len(L))
dist = lambda A, B: 3 - {1: 0, 2: 1, 4: 2, 8: 3}[len(A & B)]
cnt = {}
for B in L: cnt[dist(L[0], B)] = cnt.get(dist(L[0], B), 0) + 1
print("frames at distance 0..3 from one frame:", [cnt.get(d, 0) for d in range(4)])
def srg(n, adj):
    k = {len(adj[v]) for v in range(n)}
    lam = {len(adj[a] & adj[b]) for a, b in combinations(range(n), 2) if b in adj[a]}
    mu = {len(adj[a] & adj[b]) for a, b in combinations(range(n), 2) if b not in adj[a]}
    return n, k, lam, mu
opp = {i: {j for j in range(len(L)) if dist(L[i], L[j]) == 3} for i in range(len(L))}
print("opposite-frames graph (v,k,lambda,mu):", srg(len(L), opp))
n = 9; wt = lambda x: bin(x).count("1")
E = [x for x in range(1 << n) if wt(x) % 2 == 0]
sing = [x for x in E if x and (wt(x) // 2) % 2 == 0]
B = lambda x, y: wt(x & y) & 1
nono = {i: {j for j in range(len(sing)) if j != i and B(sing[i], sing[j])} for i in range(len(sing))}
print("singular points of O+(8,2), non-orthogonal graph (v,k,lambda,mu):", srg(len(sing), nono))
