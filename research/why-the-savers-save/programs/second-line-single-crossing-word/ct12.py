#!/usr/bin/env python3
"""an2-word: proposal 2 screen on the Coxeter-Todd lattice K12 (own code, stdlib, exact integers).

Construction (complex form): E = Z[w], w^2 + w + 1 = 0; E/2E = F_4.  K12 = { x in E^6 : x mod 2 in the hexacode }, the
hexacode being { (a, b, c, p(1), p(w), p(w^2)) : p(t) = a t^2 + b t + c } over F_4.  Checks that make it K12: minimum
norm 4, 756 minimal vectors (enumerated by brute force over all x in E^6 of norm 4, not by formula).
Screen (an-synth 4.2, as an-build's bw16.py / v-an-build's vbw16.py):
  real form: 378 lines in R^12, A = Re<x,y>/4;  complex form: 126 lines in C^6, A = <x,y>/4.
  centres t = rank, tightness (sum f f^* = (n/t) I), junk per port (other lines with A != 0), a mutually-junk clique
  (greedy, randomised restarts; a lower bound on the window k), the window of the plain label choice
  B_ST = [A_ST = 0] (rank over F_2, an upper bound on the least window for that choice), copy ratio t*k/n.
  Pass marks: copy ratio <= 0.35, window <= 12, junk per port <= 120.
usage: python3 -I ct12.py OUT.json"""
import itertools, json, random, sys
from collections import Counter
from fractions import Fraction as F

# Eisenstein integer a + b w as (a, b).  norm = a^2 - a b + b^2.
def enorm(z): a, b = z; return a * a - a * b + b * b
def emul(z, u): a, b = z; c, d = u; return (a * c - b * d, a * d + b * c - b * d)   # w^2 = -1 - w
def econj(z): a, b = z; return (a - b, -b)                                           # conj(w) = w^2 = -1 - w
def eadd(z, u): return (z[0] + u[0], z[1] + u[1])
def res2(z): return (z[0] % 2, z[1] % 2)                                             # F_4: (0,0)=0 (1,0)=1 (0,1)=w (1,1)=w^2
# F_4 arithmetic on the residues
F4 = [(0, 0), (1, 0), (0, 1), (1, 1)]
def f4add(p, q): return ((p[0] + q[0]) % 2, (p[1] + q[1]) % 2)
def f4mul(p, q): return res2(emul(p, q))
one, wf, w2f = (1, 0), (0, 1), (1, 1)
hexa = set()
for a, b, c in itertools.product(F4, repeat=3):
    def p(t): return f4add(f4add(f4mul(a, f4mul(t, t)), f4mul(b, t)), c)
    hexa.add((a, b, c, p(one), p(wf), p(w2f)))
assert len(hexa) == 64
wt = Counter(sum(1 for z in cw if z != (0, 0)) for cw in hexa)
assert wt == Counter({0: 1, 4: 45, 6: 18}), wt
# all Eisenstein integers of norm <= 4
small = [(a, b) for a in range(-3, 4) for b in range(-3, 4) if enorm((a, b)) <= 4]
bynorm = {}
for z in small: bynorm.setdefault(enorm(z), []).append(z)
assert [len(bynorm[k]) for k in (0, 1, 3, 4)] == [1, 6, 6, 6] and 2 not in bynorm
V = []
def rec(prefix, budget):
    if len(prefix) == 6:
        if budget == 0 and tuple(res2(z) for z in prefix) in hexa: V.append(tuple(prefix))
        return
    for k, zs in bynorm.items():
        if k <= budget:
            for z in zs: rec(prefix + [z], budget - k)
rec([], 4)
# minimum norm check: no nonzero lattice vector of norm 1, 2 or 3
for N in (1, 2, 3):
    W = []
    def rec2(prefix, budget):
        if len(prefix) == 6:
            if budget == 0 and tuple(res2(z) for z in prefix) in hexa: W.append(1)
            return
        for k, zs in bynorm.items():
            if k <= budget:
                for z in zs: rec2(prefix + [z], budget - k)
    rec2([], N)
    assert not W, (N, len(W))
nmin = len(V)
print('K12 (E^6, x mod 2 in the hexacode): minimum norm 4, minimal vectors', nmin)
assert nmin == 756
# hermitian product <x, y> = sum x_i conj(y_i)  (an Eisenstein integer);  real part *2 = 2a - b for a + b w
def herm(x, y):
    s = (0, 0)
    for xi, yi in zip(x, y): s = eadd(s, emul(xi, econj(yi)))
    return s
def re2(z): return 2 * z[0] - z[1]            # 2 Re(a + b w) = 2a - b
units = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]   # +-1, +-w, +-w^2  (w^2 = (-1,-1))
def scal(u, x): return tuple(emul(u, z) for z in x)
# ---------------- real form ----------------
seen, rlines = set(), []
for x in V:
    if x in seen: continue
    mx = scal((-1, 0), x); seen.add(x); seen.add(mx); rlines.append(x)
nr = len(rlines)
# real coordinates as integer vectors scaled by 2: a + b w -> (2a - b, b * sqrt3); keep the two parts separately
def rvec(x):  # returns (u, v) with real vector = (u_1/2, v_1 sqrt3/2, ...): Re-inner = (sum u u' + 3 v v')/4
    return [2 * a - b for a, b in x], [b for a, b in x]
RV = [rvec(x) for x in rlines]
def rip4(i, j):   # 4 * Re<x,y>  as an integer
    (u, v), (u2, v2) = RV[i], RV[j]
    return sum(p * q for p, q in zip(u, u2)) + 3 * sum(p * q for p, q in zip(v, v2))
assert all(rip4(i, i) == 16 for i in range(nr))     # norm 4 -> 4*4
G = [[rip4(i, j) for j in range(nr)] for i in range(nr)]
# check against the hermitian product
for i in range(0, nr, 37):
    for j in range(0, nr, 11):
        assert G[i][j] == 2 * re2(herm(rlines[i], rlines[j])), (i, j)
distr = Counter()
for i in range(nr):
    distr[tuple(sorted(Counter(abs(G[i][j]) for j in range(nr) if j != i).items()))] += 1
assert len(distr) == 1
rowdist = dict(list(distr.keys())[0])
junk_r = sum(c for k, c in rowdist.items() if k != 0)
# tightness: sum over lines of f f^T in the real 12-dim coordinates (exact with sqrt3 bookkeeping):
# coordinates (u/2, v*sqrt3/2): products u u'/4, u v' sqrt3/4, 3 v v'/4.  Need sum uu = 3 sum vv (diagonal), others 0.
Suu = [[sum(u[i] * u[j] for u, v in RV) for j in range(6)] for i in range(6)]
Svv = [[sum(v[i] * v[j] for u, v in RV) for j in range(6)] for i in range(6)]
Suv = [[sum(u[i] * v[j] for u, v in RV) for j in range(6)] for i in range(6)]
d0 = Suu[0][0]
tight_r = (all(Suu[i][j] == (d0 if i == j else 0) for i in range(6) for j in range(6)) and
           all(3 * Svv[i][j] == (d0 if i == j else 0) for i in range(6) for j in range(6)) and
           all(Suv[i][j] == 0 for i in range(6) for j in range(6)))
t_r = 12
# rank of the real Gram over Q (should be 12)
def rank_q(rows):
    M = [[F(x) for x in r] for r in rows]; rk = 0; nc = len(M[0])
    for c in range(nc):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
sample = random.Random(3).sample(range(nr), 40)
rk_r = rank_q([[G[i][j] for j in range(nr)] for i in sample])
def rank_f2(rows):
    basis = {}
    for r in rows:
        while r:
            h = r.bit_length() - 1
            if h in basis: r ^= basis[h]
            else: basis[h] = r; break
    return len(basis)
def clique(adj, n, trials, seed):
    rnd = random.Random(seed); best = []
    for _ in range(trials):
        order = list(range(n)); rnd.shuffle(order); cl = []; cand = set(range(n))
        for i in order:
            if i in cand: cl.append(i); cand &= adj[i]
        if len(cl) > len(best): best = cl
    return best
adj_r = {i: {j for j in range(nr) if j != i and G[i][j] != 0} for i in range(nr)}
cl_r = clique(adj_r, nr, 400, 1)
assert all(G[a][b] != 0 for a, b in itertools.combinations(cl_r, 2))
B_r = [sum(1 << j for j in range(nr) if j == i or G[i][j] == 0) for i in range(nr)]
k_r = rank_f2(B_r)
real = dict(lines=nr, centres_t=t_r, gram_rank_on_40_sampled_rows=rk_r, n_over_t=F(nr, t_r).__str__(), tight=tight_r,
            abs_4Re_ip_from_each_line=rowdist, normalised='|A| = |4Re<x,y>|/16: 4 -> 1/4, 8 -> 1/2',
            junk_per_port=junk_r, junk_clique_found=len(cl_r), window_plain_labels_rankF2=k_r,
            copy_ratio_lower=t_r * len(cl_r) / nr, copy_ratio_plain=t_r * k_r / nr,
            one_plus_sum_A2=str(1 + sum(F(c) * F(k, 16) ** 2 for k, c in rowdist.items() if k != 0)))
real['passes'] = dict(copy_ratio=real['copy_ratio_lower'] <= 0.35, window=len(cl_r) <= 12, junk=junk_r <= 120)
print('REAL form:', json.dumps(real, indent=1))
# ---------------- complex form ----------------
seen, clines = set(), []
for x in V:
    if x in seen: continue
    for u in units: seen.add(scal(u, x))
    clines.append(x)
nc = len(clines)
H = [[herm(clines[i], clines[j]) for j in range(nc)] for i in range(nc)]
cdist = Counter()
for i in range(nc):
    cdist[tuple(sorted(Counter(enorm(H[i][j]) for j in range(nc) if j != i).items()))] += 1
assert len(cdist) == 1
crow = dict(list(cdist.keys())[0])
junk_c = sum(c for k, c in crow.items() if k != 0)
adj_c = {i: {j for j in range(nc) if j != i and H[i][j] != (0, 0)} for i in range(nc)}
cl_c = clique(adj_c, nc, 400, 2)
B_c = [sum(1 << j for j in range(nc) if j == i or H[i][j] == (0, 0)) for i in range(nc)]
k_c = rank_f2(B_c)
t_c = 6
comp = dict(lines=nc, centres_t=t_c, n_over_t=str(F(nc, t_c)), abs2_ip_from_each_line=crow,
            normalised='|A|^2 = |<x,y>|^2/16: 4 -> 1/4, 3 -> 3/16, 1 -> 1/16',
            junk_per_port=junk_c, junk_clique_found=len(cl_c), window_plain_labels_rankF2=k_c,
            copy_ratio_lower=t_c * len(cl_c) / nc, copy_ratio_plain=t_c * k_c / nc,
            one_plus_sum_A2=str(1 + sum(F(c) * F(k, 16) for k, c in crow.items() if k != 0)),
            phases_note='entries of A are Eisenstein (powers of w = exp(2 pi i/3)), not powers of i')
comp['passes'] = dict(copy_ratio=comp['copy_ratio_lower'] <= 0.35, window=len(cl_c) <= 12, junk=junk_c <= 120)
print('COMPLEX form:', json.dumps(comp, indent=1))
json.dump(dict(minimal_vectors=nmin, real=real, complex=comp, real_clique=[rlines[i] for i in cl_r]),
          open(sys.argv[1], 'w'), indent=1)
