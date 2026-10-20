#!/usr/bin/env python3
"""v-an2-word: independent Coxeter-Todd screen (own code, stdlib, exact integers).

Construction DIFFERENT from an2-word's (which used x mod 2 in the hexacode):
  E = Z[w], theta = w - w^2 = sqrt(-3).  K12' = { x in E^6 : x_1 = x_2 = ... = x_6 (mod theta), sum x_i = 0 (mod 3) }.
  Expected: minimum norm 6, 756 minimal vectors (270 of shape (theta u, -theta u', 0^4), 486 of shape (u_1..u_6)).
  Minimal vectors are enumerated by brute force over all x in E^6 of norm <= 6 (no formula), and norms 1..5 are
  checked to be empty.
Screen quantities (as an-build bw16.py / an2-word ct12.py): lines (real: up to sign; complex: up to the 6 units),
  inner-product distribution from each line, junk per port (other lines with nonzero product), tightness of the frame
  operator, a mutually-junk clique (greedy with restarts, plus a structured clique), copy ratio t*k/n.
usage: python3 -I vct12.py OUT.json"""
import itertools, json, random, sys
from collections import Counter

def enorm(z): a, b = z; return a * a - a * b + b * b
def emul(z, u): a, b = z; c, d = u; return (a * c - b * d, a * d + b * c - b * d)   # w^2 = -1 - w
def econj(z): a, b = z; return (a - b, -b)
def eadd(z, u): return (z[0] + u[0], z[1] + u[1])
def esub(z, u): return (z[0] - u[0], z[1] - u[1])
THETA = (1, 2)            # w - w^2 = w - (-1 - w) = 1 + 2w ; norm 1 - 2 + 4 = 3
assert enorm(THETA) == 3
def divisible(z, d):
    # z / d in E ?  z * conj(d) / N(d)
    p = emul(z, econj(d)); n = enorm(d)
    return p[0] % n == 0 and p[1] % n == 0
def cong(z, u, d): return divisible(esub(z, u), d)
small = [(a, b) for a in range(-4, 5) for b in range(-4, 5) if enorm((a, b)) <= 6]
bynorm = {}
for z in small: bynorm.setdefault(enorm(z), []).append(z)
print('Eisenstein norms <= 6 and counts:', {k: len(v) for k, v in sorted(bynorm.items())})
def inlat(x):
    s = (0, 0)
    for z in x: s = eadd(s, z)
    if not divisible(s, (3, 0)): return False
    return all(cong(x[i], x[0], THETA) for i in range(1, 6))
def enum(N):
    out = []
    def rec(prefix, budget):
        if len(prefix) == 6:
            if budget == 0 and inlat(prefix): out.append(tuple(prefix))
            return
        for k, zs in bynorm.items():
            if k <= budget:
                for z in zs: rec(prefix + [z], budget - k)
    rec([], N)
    return out
for N in range(1, 6):
    assert not enum(N), N
V = enum(6)
print('minimal norm 6, minimal vectors:', len(V))
assert len(V) == 756
shape = Counter(tuple(sorted(enorm(z) for z in x)) for x in V)
print('shapes (sorted coordinate norms):', dict(shape))
# lattice determinant check (sanity): index computations skipped; instead check closure of V under the 6 units
units = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]
Vs = set(V)
assert all(tuple(emul(u, z) for z in x) in Vs for x in V for u in units)
def herm(x, y):
    s = (0, 0)
    for xi, yi in zip(x, y): s = eadd(s, emul(xi, econj(yi)))
    return s
def re2(z): return 2 * z[0] - z[1]            # 2 Re(a + b w)
# real form
seen, R = set(), []
for x in V:
    if x in seen: continue
    seen.add(x); seen.add(tuple(emul((-1, 0), z) for z in x)); R.append(x)
nr = len(R)
G = [[re2(herm(R[i], R[j])) for j in range(nr)] for i in range(nr)]      # 2 Re<x,y>, norm -> 12
assert all(G[i][i] == 12 for i in range(nr))
rows = Counter()
for i in range(nr):
    rows[tuple(sorted(Counter(abs(G[i][j]) for j in range(nr) if j != i).items()))] += 1
assert len(rows) == 1, rows
rrow = dict(list(rows)[0])
junk_r = sum(c for k, c in rrow.items() if k)
# tightness in real coordinates: z = a + b w -> (a - b/2, b sqrt3/2).  Frame operator sum v v^T must be c I.
# Use scaled coordinates p = 2a - b, q = b (real vector = (p/2, q sqrt3/2)): need sum p p' = 3 sum q q' on the diagonal,
# zero off-diagonal, and sum p q' = 0.
P = [[2 * a - b for a, b in x] for x in R]; Q = [[b for a, b in x] for x in R]
Spp = [[sum(p[i] * p[j] for p in P) for j in range(6)] for i in range(6)]
Sqq = [[sum(q[i] * q[j] for q in Q) for j in range(6)] for i in range(6)]
Spq = [[sum(p[i] * q[j] for p, q in zip(P, Q)) for j in range(6)] for i in range(6)]
d0 = Spp[0][0]
tight_r = (all(Spp[i][j] == (d0 if i == j else 0) for i in range(6) for j in range(6)) and
           all(3 * Sqq[i][j] == (d0 if i == j else 0) for i in range(6) for j in range(6)) and
           all(Spq[i][j] == 0 for i in range(6) for j in range(6)))
def greedy_clique(adj, n, trials, seed):
    rnd = random.Random(seed); best = []
    for _ in range(trials):
        order = list(range(n)); rnd.shuffle(order); cl = []; cand = set(range(n))
        for i in order:
            if i in cand: cl.append(i); cand &= adj[i]
        if len(cl) > len(best): best = cl
    return best
adj_r = {i: {j for j in range(nr) if j != i and G[i][j] != 0} for i in range(nr)}
cl_r = greedy_clique(adj_r, nr, 2000, 11)
assert all(G[a][b] != 0 for a, b in itertools.combinations(cl_r, 2))
# complex form
seen, C = set(), []
for x in V:
    if x in seen: continue
    for u in units: seen.add(tuple(emul(u, z) for z in x))
    C.append(x)
nc = len(C)
H = [[enorm(herm(C[i], C[j])) for j in range(nc)] for i in range(nc)]   # |<x,y>|^2, norm -> 36
assert all(H[i][i] == 36 for i in range(nc))
crows = Counter()
for i in range(nc):
    crows[tuple(sorted(Counter(H[i][j] for j in range(nc) if j != i).items()))] += 1
assert len(crows) == 1
crow = dict(list(crows)[0])
junk_c = sum(c for k, c in crow.items() if k)
adj_c = {i: {j for j in range(nc) if j != i and H[i][j] != 0} for i in range(nc)}
cl_c = greedy_clique(adj_c, nc, 4000, 12)
assert all(H[a][b] != 0 for a, b in itertools.combinations(cl_c, 2))
# complex tightness: sum x x^* = c I (hermitian, exact over E)
Sc = [[(0, 0)] * 6 for _ in range(6)]
for x in C:
    for i in range(6):
        for j in range(6): Sc[i][j] = eadd(Sc[i][j], emul(x[i], econj(x[j])))
tight_c = all(Sc[i][j] == (Sc[0][0] if i == j else (0, 0)) for i in range(6) for j in range(6))
# products of a complex line pair: which phases occur (Eisenstein units times what?)
phases = Counter()
for j in range(1, nc):
    z = herm(C[0], C[j])
    if z != (0, 0): phases[z] += 1
out = dict(construction='x_i all congruent mod theta, sum = 0 mod 3, in E^6; min norm 6',
           minimal_vectors=len(V), shapes={str(k): v for k, v in shape.items()},
           real=dict(lines=nr, t=12, n_over_t=nr / 12, tight=tight_r,
                     abs_2Re_ip_from_each_line={str(k): v for k, v in sorted(rrow.items())},
                     normalised='|2Re<x,y>|/12: 3 -> 1/4, 6 -> 1/2', junk_per_port=junk_r,
                     junk_clique_found=len(cl_r), copy_ratio_lower=12 * len(cl_r) / nr,
                     passes=dict(copy_ratio=12 * len(cl_r) / nr <= 0.35, window=len(cl_r) <= 12, junk=junk_r <= 120)),
           complex=dict(lines=nc, t=6, n_over_t=nc / 6, tight=tight_c,
                        abs2_ip_from_each_line={str(k): v for k, v in sorted(crow.items())},
                        normalised='|<x,y>|^2/36: 9 -> 1/4', junk_per_port=junk_c,
                        junk_clique_found=len(cl_c), copy_ratio_lower=6 * len(cl_c) / nc,
                        passes=dict(copy_ratio=6 * len(cl_c) / nc <= 0.35, window=len(cl_c) <= 12, junk=junk_c <= 120),
                        products_from_line0={str(k): v for k, v in phases.items()}))
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
