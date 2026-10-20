"""cover_exact.py (v-g-9hyp-model) -- exact re-check of d-bitspan's F1, F2, F3 with certificates that do not rely on
float ranks.  F1: exact rational rank of r(q) for all 511 q != 0 in W0.  F2: min |r(q)|.  F3: for every q != 0 and
every root-orthogonal 6-flat F (perp of a pair of lines: 3780 D6 + 1120 E6), rank_Q(r(q) minus F) >= 7, certified by
rank mod p of the Gram matrix (rank_p(M^T M) <= rank_Q(M)).  usage: python -I cover_exact.py <s2|word>"""
import sys, time, collections; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from fractions import Fraction
import numpy as np
import geo9
name = sys.argv[1]
P, live, W0b = geo9.unit(name); v = len(live); LA = geo9.LA; G = geo9.GR
allW0 = [x for x in geo9.f2_span_all(W0b) if x]
assert len(allW0) == 511 and len(set(allW0)) == 511
def qrank(rows):
    M = [[Fraction(int(a)) for a in r] for r in rows]; rk = 0
    for c in range(8):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
supp = [[t for t in range(v) if (x >> t) & 1] for x in allW0]
F1 = collections.Counter(qrank([LA[live[t]] for t in s]) for s in supp)
print(name, 'F1 exact Q-rank of r(q), q != 0:', sorted(F1.items()))
print(name, 'F2 min |r(q)| =', min(len(s) for s in supp), ' sizes>=44 needed for one part >= 22')
# flats
fl = {}
for l in range(120):
    for m in range(l + 1, 120):
        key = frozenset(T for T in range(120) if G[T, l] == 0 and G[T, m] == 0)
        fl.setdefault(key, []).append((l, m))
print(name, 'distinct perp-of-pair flats', len(fl), 'line counts', sorted(collections.Counter(len(k) for k in fl).items()),
      'exact ranks', sorted(collections.Counter(qrank([LA[T] for T in k]) for k in list(fl)[:200]).items()), '(first 200)')
Fm = np.array([[1 if live[t] in k else 0 for t in range(v)] for k in fl], dtype=np.int64)   # nF x v
R = np.array([[1 if (x >> t) & 1 else 0 for t in range(v)] for x in allW0], dtype=np.int64)  # 511 x v
LL = np.array([np.outer(LA[T], LA[T]).reshape(64) for T in live], dtype=np.int64)       # v x 64
p = 2147483647
def inv_mod(a):
    r = np.ones_like(a); b = a.copy(); e = p - 2
    while e:
        if e & 1: r = (r * b) % p
        b = (b * b) % p; e >>= 1
    return r
def batch_rank(A):
    A = A % p; N = A.shape[0]; used = np.zeros((N, 8), bool); rk = np.zeros(N, np.int64); ar = np.arange(N)
    for c in range(8):
        cand = (A[:, :, c] != 0) & ~used
        has = cand.any(1); piv = np.argmax(cand, 1)
        pr = A[ar, piv]                                   # N x 8
        iv = inv_mod(np.where(has, pr[:, c], 1))
        f = (A[:, :, c] * iv[:, None]) % p
        f[ar, piv] = 0; f[~has] = 0
        A = (A - (f[:, :, None] * pr[:, None, :]) % p) % p
        used[ar[has], piv[has]] = True; rk += has
    return rk
t0 = time.time(); hist = collections.Counter(); low = []
for qi in range(511):
    rest = R[qi][None, :] * (1 - Fm)                      # lines of r(q) outside F
    Gm = (rest @ LL).reshape(-1, 8, 8)
    rk = batch_rank(Gm)
    m = int(rk.min()); hist[m] += 1
    if m <= 6: low.append((qi, m, int((rk <= 6).sum())))
print(name, 'F3 min over F of rank_p(Gram(r(q) minus F)), histogram over q:', sorted(hist.items()), ' q with <= 6:', low[:10],
      ' %.1f s' % (time.time() - t0))
# sanity of the mod-p rank routine on random subsets vs exact
rng = np.random.default_rng(1); bad = 0
for _ in range(300):
    k = int(rng.integers(1, 12)); S = rng.choice(v, size=k, replace=False)
    g = sum(LL[t] for t in S).reshape(1, 8, 8)
    if int(batch_rank(g)[0]) != qrank([LA[live[t]] for t in S]): bad += 1
print(name, 'mod-p Gram rank vs exact Q rank mismatches on 300 random subsets:', bad)
