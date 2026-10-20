# l3n.py: Corollary L3 (block 3, Singer Z_7) and Proposition N orbit ratios (blocks 3, 4), own code.
import itertools
from fractions import Fraction as Fr
def setup(m):
    P = [(i, j) for i in range(m) for j in range(i, m)]
    D = list(range(1, 1 << m))
    v = lambda x: [(x >> i) & 1 for i in range(m)]
    q = lambda Q, x: sum(v(x)[i] * v(x)[j] for b, (i, j) in enumerate(P) if (Q >> b) & 1) % 2
    pr = lambda Q, S: sum(S[i][j] for b, (i, j) in enumerate(P) if (Q >> b) & 1) % 2
    mul = lambda A, B: tuple(tuple(sum(A[i][k] * B[k][j] for k in range(m)) % 2 for j in range(m)) for i in range(m))
    tp = lambda A: tuple(tuple(A[j][i] for j in range(m)) for i in range(m))
    return P, D, v, q, pr, mul, tp
# --- L3
m = 3; P, D, v, q, pr, mul, tp = setup(m)
def gm(a, b):
    r = 0
    for i in range(3):
        if (b >> i) & 1: r ^= a << i
    for d in (4, 3):
        if (r >> d) & 1: r ^= 0b1011 << (d - 3)
    return r
cols = [gm(2, 1 << j) for j in range(3)]
A = tuple(tuple((cols[j] >> i) & 1 for j in range(3)) for i in range(3))
ID = tuple(tuple(int(i == j) for j in range(3)) for i in range(3))
H = [ID]; B = A
while B != ID: H.append(B); B = mul(A, B)
T = [mul(h, tp(h)) for h in H]; print('block 3 Singer |H| =', len(H), 'distinct targets', len(set(T)))
def form(s):
    Q = 0
    for t in s.replace(' ', '').split('+'):
        w = [int(c) - 1 for c in t.split('x') if c]
        i, j = (w[0], w[0]) if len(w) == 1 else (min(w), max(w)); Q |= 1 << P.index((i, j))
    return Q
seven = ['x1 + x1x3', 'x1x2 + x1x3', 'x2 + x1x2', 'x1 + x1x2 + x1x3 + x2x3', 'x2 + x1x2 + x2x3', 'x3 + x2x3', 'x3 + x1x3 + x2x3']
cover = {x: 0 for x in D}
for s in seven:
    Q = form(s); S = [x for x in D if q(Q, x)]
    for x in S: cover[x] += 1
    print('%-28s supp %s odd at %d of 7 targets' % (s, S, sum(pr(Q, t) for t in T)))
print('each direction covered', set(cover.values()), 'times')
# --- Proposition N orbit ratios
for m in (3, 4):
    P, D, v, q, pr, mul, tp = setup(m)
    I = tuple(tuple(int(i == j) for j in range(m)) for i in range(m))
    def det(A):
        M_ = [list(r) for r in A]
        for c in range(m):
            p = next((r for r in range(c, m) if M_[r][c]), None)
            if p is None: return 0
            M_[c], M_[p] = M_[p], M_[c]
            for r in range(m):
                if r != c and M_[r][c]: M_[r] = [(a + b) % 2 for a, b in zip(M_[r], M_[c])]
        return 1
    GL = [A for A in (tuple(tuple((r[i] >> j) & 1 for j in range(m)) for i in range(m)) for r in itertools.product(range(1 << m), repeat=m)) if det(A)]
    act = lambda A, x: sum((sum(A[i][j] * v(x)[j] for j in range(m)) % 2) << i for i in range(m))
    NF = 1 << len(P); supp = {Q: frozenset(x for x in D if q(Q, x)) for Q in range(1, NF)}
    key = {s: Q for Q, s in supp.items()}
    seen = set(); out = []
    for Q in range(1, NF):
        if Q in seen: continue
        orb = {key[frozenset(act(g, x) for x in supp[Q])] for g in GL}; seen |= orb
        a = Fr(sum(pr(Q2, I) for Q2 in orb), len(orb)); s = len(supp[Q])
        out.append((len(orb), s, a, float((2 ** m - 1) * a / s)))
    print('block', m, 'orbits (size, |supp|, a_O at I, (2^m-1) a/s):', [(o[0], o[1], str(o[2]), round(o[3], 3)) for o in out], 'max', max(o[3] for o in out))
