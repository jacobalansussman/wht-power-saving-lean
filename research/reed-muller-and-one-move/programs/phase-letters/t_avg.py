#!/usr/bin/env python3
"""t_avg.py (rm2-phase): the Averaging Lemma written directly in Lemma N's block C[Sp(4,2)] (block 2), as a code check.
Enumerate Sp(4,2) (720 elements) and P = Stab(L0) (48).  Take a random family Y_U (U = the 15 frames) of subspaces of
Q^4 with unequal dimensions and no symmetry.  In V = Map(Sp, Q^4) (dimension 2880, (g v)_gam = v_(g^-1 gam)) put
Yhat = {v : v_gam in Y_(gam^-1 L0)}.  Check: (a) p Yhat = Yhat for all p in P (by recomputing the frame labels);
(b) f(g) = dim(Yhat + g Yhat) - dim Yhat equals |Sp| * abar_d(gL0, L0) for EVERY g (720 checks);
(c) so moves t_z have f = 720 abar_1 and K has f = 720 abar_2.  Ranks exact (sympy-free integer Bareiss)."""
import random, itertools
from spg import Sp
from fractions import Fraction as Fr

def rank_int(rows):
    M = [list(r) for r in rows if any(r)]
    if not M: return 0
    n = len(M[0]); r = 0; prev = 1
    for c in range(n):
        p = next((k for k in range(r, len(M)) if M[k][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for k in range(r + 1, len(M)):
            M[k] = [(M[r][c] * M[k][j] - M[k][c] * M[r][j]) // prev for j in range(n)]
        prev = M[r][c]; r += 1
        if r == len(M): break
    return r

sp = Sp(2); rng = random.Random(7)
# enumerate the group
gens = [sp.t(z) for z in (1, 2, 3)] + [sp.yC([1, 0]), sp.yC([0, 2]), sp.yC([2, 1]), sp.lev([1 | 2, 2]), sp.lev([2, 1])]
G = {sp.ident()}; st = [sp.ident()]
while st:
    g = st.pop()
    for h in gens:
        k = sp.mul(h, g)
        if k not in G: G.add(k); st.append(k)
G = sorted(G); L0 = sp.L0()
P = [g for g in G if sp.act(g, L0) == L0]
frames = sorted({sp.act(g, L0) for g in G}, key=lambda U: sorted(U))
print("|Sp(4,2)| =", len(G), " |P| =", len(P), " frames", len(frames))
# random family
Y = {}
for U in frames:
    d = rng.choice([1, 2, 2, 3])
    Y[U] = [[rng.randint(-2, 2) for _ in range(4)] for _ in range(d)]
dimY = {U: rank_int(Y[U]) for U in frames}
def delta(A, B):   # dim Y_A - dim(Y_A cap Y_B) = dim(Y_A + Y_B) - dim Y_B
    return rank_int(Y[A] + Y[B]) - dimY[B]
abar = {}
for A in frames:
    for B in frames:
        d = sp.dist(A, B)
        if d == 0: continue
        abar.setdefault(d, []).append(delta(A, B))
abar = {d: Fr(sum(v), len(v)) for d, v in abar.items()}
print("abar_1 =", abar[1], " abar_2 =", abar[2], " (pairs:", {d: None for d in abar}, ")")
# (a) keptness: the label of position gam in p.Yhat is (p^-1 gam)^-1 L0 = gam^-1 p L0 = gam^-1 L0
inv = {g: sp.inv(g) for g in G}
oka = all(sp.act(inv[sp.mul(inv[p], gam)], L0) == sp.act(inv[gam], L0) for p in P for gam in G)
# (b) f(g) position by position: position gam holds Y_(gam^-1 L0) in Yhat and Y_(gam^-1 g L0) in g.Yhat
okb = True; seen = {}
for g in G:
    gL0 = sp.act(g, L0); d = sp.dist(gL0, L0)
    f = 0
    for gam in G:
        A = sp.act(inv[gam], gL0); Bf = sp.act(inv[gam], L0)
        f += rank_int(Y[Bf] + Y[A]) - dimY[Bf]
    target = 0 if d == 0 else len(G) * abar[d]
    okb &= (f == target); seen[d] = seen.get(d, 0) + 1
print("(a) every p in P keeps Yhat:", oka)
print("(b) f(g) = |Sp| * abar_d for all", len(G), "elements (distance counts", dict(sorted(seen.items())), "):", okb)
t = sp.t(1); K = sp.K()
print("(c) f(t_e1) = 720 abar_1 =", 720 * abar[1], "; f(K) = 720 abar_2 =", 720 * abar[2],
      "; Theorem R in C[Sp(4,2)]: moves >= W * abar_2 / abar_1 =", abar[2] / abar[1], "W")
