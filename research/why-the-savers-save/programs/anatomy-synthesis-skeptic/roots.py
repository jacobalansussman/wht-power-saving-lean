"""v-an-synth: independent check of the E8 / E7 design numbers (own code, stdlib, written from scratch;
does not import or read an-synth's calc).  Run: python3 -I calc/roots.py"""
import itertools, json, sys
from fractions import Fraction as F

def e8_roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * 8; v[i] = 2 * si; v[j] = 2 * sj; R.append(tuple(v))   # doubled coords
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.append(tuple(signs))   # doubled: (+-1/2)^8 -> (+-1)^8
    return R

def ip(a, b):  # true inner product = sum(a*b)/4 for doubled coords
    s = sum(x * y for x, y in zip(a, b)); assert s % 4 == 0; return s // 4

def lines(R):
    L = []; seen = set()
    for r in R:
        key = max(r, tuple(-x for x in r))
        if key not in seen:
            seen.add(key); L.append(key)
    return L

def rank_f2(rows):
    rows = [int(''.join('1' if b else '0' for b in r), 2) for r in rows]
    basis = []
    for x in rows:
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    return len(basis)

def rank_q(M):
    M = [[F(x) for x in row] for row in M]; n = len(M); m = len(M[0]); r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

def max_clique(adj, n):
    best = []
    def bk(R, P, X):
        nonlocal best
        if not P and not X:
            if len(R) > len(best): best = R[:]
            return
        if len(R) + len(P) <= len(best): return
        u = max(P | X, key=lambda w: len(adj[w] & P))
        for v in list(P - adj[u]):
            bk(R + [v], P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
    bk([], set(range(n)), set())
    return best

def design(name, L):
    n = len(L)
    G = [[ip(a, b) for b in L] for a in L]           # norm 2 roots: entries 2, +-1, 0
    t = rank_q(G)
    junk = [sum(1 for j in range(n) if j != i and G[i][j] != 0) for i in range(n)]
    B = [[(1 + G[i][j]) % 2 for j in range(n)] for i in range(n)]
    k = rank_f2(B)
    # tightness A^2 = (n/t) A with A = G/2
    A = [[F(x, 2) for x in row] for row in G]
    tight = all(sum(A[i][l] * A[l][j] for l in range(n)) == F(n, t) * A[i][j] for i in range(min(n, 6)) for j in range(n))
    junk_adj = [set(j for j in range(n) if j != i and G[i][j] != 0) for i in range(n)]
    orth_adj = [set(j for j in range(n) if j != i and G[i][j] == 0) for i in range(n)]
    om = max_clique(junk_adj, n); oo = max_clique(orth_adj, n)
    out = dict(design=name, n=n, t=t, junk_min=min(junk), junk_max=max(junk), n_over_t=n / t,
               one_plus_s4=1 + sum(junk) / n / 4, tight_rows_checked=tight, window_rank_f2_JplusG=k,
               junk_clique=len(om), orth_clique=len(oo), copy_ratio=t * k / n,
               copy_ratio_floor_with_clique_window=t * len(om) / n)
    print(out); return out

R = e8_roots(); assert len(R) == 240
L8 = lines(R)
res = [design('E8', L8)]
r0 = L8[0]
L7 = [l for l in L8 if ip(l, r0) == 0]
res.append(design('E7 (E8 lines orthogonal to one root)', L7))
# P112: E8 minus an orthogonal frame of 8 root lines e1+-e2, e3+-e4, e5+-e6, e7+-e8 (doubled coords)
frame = set()
for a in (0, 2, 4, 6):
    for s in (1, -1):
        v = [0] * 8; v[a] = 2; v[a + 1] = 2 * s; frame.add(max(tuple(v), tuple(-x for x in v)))
P112 = [l for l in L8 if l not in frame]
assert len(P112) == 112
res.append(design('P112 = E8 minus an orthogonal frame', P112))
json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else 'roots.json', 'w'), indent=1)
