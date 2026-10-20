#!/usr/bin/env python3
"""v_star.py (v-e3-slice6): exact tests of the star lemma over GF(p) (the lemma claims to hold over any field).
Statement tested: E nonsingular on X; X = Psi + Q_1 + ... + Q_g (disjoint); R, C in Psi, |R|=|C|=k, E[R,C] nonsingular;
sigma_j = rank E[Psi u Q_j, Psi u Q_j]; a = rank of the Schur complement T = E/E[R,C] on (Psi\\R) x (Psi\\C). Then
   2(|Psi| - k) <= sum_j (sigma_j - k) - (g - 2) a.
Instances: (i) AG(4,2) flats (Psi a d-flat, pieces the (d+1)-flats through it), (ii) arbitrary random partitions,
matrices random over GF(2), GF(3), GF(10007), plus low-rank-structured ones and the bent matrix mod a big prime.
"""
import random, itertools, sys
def rank_mod(M, p):
    M = [row[:] for row in M]; r = 0
    rows = len(M); cols = len(M[0]) if rows else 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c] % p: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r
def inv_mod(M, p):
    n = len(M); A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if A[i][c] % p)
        A[c], A[piv] = A[piv], A[c]
        iv = pow(A[c][c], p - 2, p); A[c] = [(x * iv) % p for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] % p:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]
def sub(E, R, C): return [[E[r][c] for c in C] for r in R]
def schur(E, X, R, C, p):
    Rr = [x for x in X if x not in R]; Cc = [x for x in X if x not in C]
    Ai = inv_mod(sub(E, R, C), p)
    B = sub(E, Rr, C); Cm = sub(E, R, Cc); D = sub(E, Rr, Cc)
    # T = D - B Ai Cm
    BA = [[sum(B[i][t] * Ai[t][j] for t in range(len(R))) % p for j in range(len(R))] for i in range(len(Rr))]
    T = {}
    for i, r in enumerate(Rr):
        for j, c in enumerate(Cc):
            T[(r, c)] = (D[i][j] - sum(BA[i][t] * Cm[t][j] for t in range(len(R)))) % p
    return T
n = 4; X = list(range(16))
def pc(x): return bin(x).count('1')
def flats(d):
    out = set()
    for basis in itertools.combinations(range(1, 16), d):
        span = {0}
        for b in basis: span |= {s ^ b for s in span}
        if len(span) != 2 ** d: continue
        for a in range(16): out.add(frozenset(a ^ s for s in span))
    return sorted(out, key=sorted)
FL = {d: flats(d) for d in range(0, 4)}
def pieces_flats(Psi, d):
    return [sorted(G - Psi) for G in FL[d + 1] if Psi <= G]
def random_matrix(p, kind):
    while True:
        if kind == 'rand':
            E = [[random.randrange(p) for _ in X] for _ in X]
        elif kind == 'lowrank':  # sum of a few rank-1 pieces supported on random flats, plus a permutation
            E = [[0] * 16 for _ in X]
            perm = list(range(16)); random.shuffle(perm)
            for i in X: E[i][perm[i]] = random.randrange(1, p)
            for _ in range(random.randrange(1, 6)):
                d = random.choice([1, 2, 3]); G = random.choice(FL[d]); H = random.choice(FL[d])
                u = {x: random.randrange(p) for x in G}; v = {y: random.randrange(p) for y in H}
                for x in G:
                    for y in H: E[x][y] = (E[x][y] + u[x] * v[y]) % p
        if rank_mod(E, p) == 16: return E
def bent(p):
    phi = lambda x: ((x & 1) * ((x >> 1) & 1) + ((x >> 2) & 1) * ((x >> 3) & 1)) & 1
    return [[(1 if phi(x ^ y) == 0 else p - 1) for y in X] for x in X]
def test(E, Psi, pieces, p, k):
    Psi = sorted(Psi)
    # random R, C of size k inside Psi with E[R,C] nonsingular
    for _ in range(30):
        R = random.sample(Psi, k); C = random.sample(Psi, k)
        if k == 0 or rank_mod(sub(E, R, C), p) == k: break
    else: return None
    sig = [rank_mod(sub(E, Psi + Qj, Psi + Qj), p) for Qj in pieces]
    if k:
        T = schur(E, X, R, C, p)
        A = [[T[(r, c)] for c in Psi if c not in C] for r in Psi if r not in R]
        a = rank_mod(A, p) if A and A[0] else 0
    else:
        a = rank_mod(sub(E, Psi, Psi), p)
    g = len(pieces)
    lhs = 2 * (len(Psi) - k); rhs = sum(s - k for s in sig) - (g - 2) * a
    return lhs, rhs
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
viol = 0; tight = 0; tot = 0
for trial in range(1500):
    p = random.choice([2, 3, 10007]); kind = random.choice(['rand', 'lowrank', 'lowrank'])
    E = random_matrix(p, kind) if trial % 10 else bent(1000003 if trial % 20 else 3)
    if trial % 10 == 0: p = 1000003 if trial % 20 else 3
    if trial % 2:
        d = random.choice([0, 1, 2]); Psi = set(random.choice(FL[d])); pieces = pieces_flats(frozenset(Psi), d)
    else:   # arbitrary partition
        perm = X[:]; random.shuffle(perm); sz = random.randrange(1, 8); Psi = set(perm[:sz]); rest = perm[sz:]
        g = random.randrange(1, 6); cuts = sorted(random.sample(range(1, len(rest)), min(g - 1, len(rest) - 1))) if len(rest) > 1 else []
        pieces = [rest[i:j] for i, j in zip([0] + cuts, cuts + [len(rest)])]
    rk = rank_mod(sub(E, sorted(Psi), sorted(Psi)), p)
    k = random.randrange(0, rk + 1)
    res = test(E, Psi, pieces, p, k)
    if res is None: continue
    tot += 1; lhs, rhs = res
    if lhs > rhs: viol += 1; print("VIOLATION", p, kind, sorted(Psi), pieces, k, lhs, rhs)
    if lhs == rhs: tight += 1
print("star lemma tests:", tot, "violations:", viol, "tight:", tight)
# bent matrix exact facts over a big prime (the theorem's 'L(4,4) feasible' and the plane ranks)
p = 1000003; E = bent(p)
print("bent mod p: rank", rank_mod(E, p), "; hyperplane ranks", sorted(set(rank_mod(sub(E, sorted(H), sorted(H)), p) for H in FL[3])),
      "; plane ranks", sorted(set(rank_mod(sub(E, sorted(P), sorted(P)), p) for P in FL[2])))
# --- diagnostics: distribution of (|Psi|, g, k, lhs, rhs) on the bent matrix with flats, and on random low-rank
from collections import Counter
cnt = Counter(); random.seed(7)
for t in range(400):
    p = 1000003; E = bent(p); d = random.choice([0, 1, 2]); Psi = frozenset(random.choice(FL[d]))
    pieces = pieces_flats(Psi, d); rk = rank_mod(sub(E, sorted(Psi), sorted(Psi)), p); k = random.randrange(0, rk + 1)
    res = test(E, set(Psi), pieces, p, k)
    if res: cnt[(len(Psi), len(pieces), k) + res] += 1
for key in sorted(cnt): print("bent  |Psi|,g,k,lhs,rhs =", key, "x", cnt[key])
cnt = Counter()
for t in range(400):
    p = 2; E = random_matrix(p, 'lowrank'); d = random.choice([0, 1, 2]); Psi = frozenset(random.choice(FL[d]))
    pieces = pieces_flats(Psi, d); rk = rank_mod(sub(E, sorted(Psi), sorted(Psi)), p); k = random.randrange(0, rk + 1)
    res = test(E, set(Psi), pieces, p, k)
    if res: cnt[(len(Psi), len(pieces), k) + res] += 1
print("GF(2) low-rank: min slack rhs-lhs by (|Psi|,g,k):", sorted({(a, b, c): None for (a, b, c, l, r) in cnt}.keys())[:3],
      "; min rhs-lhs =", min(r - l for (a, b, c, l, r) in cnt), "; cases", sum(cnt.values()))
