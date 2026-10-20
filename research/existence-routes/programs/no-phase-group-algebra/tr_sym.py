#!/usr/bin/env python3
"""tr_sym.py (e3-rmbuild): Theorem R certificates for the no-phase game via Theorem A.
V = C[Sym_m] (Sym_m acting by translation, GL(m,2) by congruence S -> h S h^T); Y_0 = C[U] for U a union of
GL-orbits on Sym_m (GL-invariant, so every free letter keeps Y_0).  f(g) = dim(Y_0 + g Y_0) - dim Y_0 = |(U+g)\\U|.
Theorem R: M * f(zz^T) >= W * f(I)  (all moves have the same f since GL is transitive on rank-one matrices).
Prints the best ratio f(I)/f(zz^T) over all unions of orbits (and the orbits)."""
import itertools, sys
m = int(sys.argv[1]) if len(sys.argv) > 1 else 4
idx = [(j, k) for j in range(m) for k in range(j, m)]
def enc(S):          # S as m x m 0/1 list -> int
    v = 0
    for b, (j, k) in enumerate(idx):
        if S[j][k]: v |= 1 << b
    return v
def dec(v):
    S = [[0] * m for _ in range(m)]
    for b, (j, k) in enumerate(idx):
        if (v >> b) & 1: S[j][k] = S[k][j] = 1
    return S
nS = 1 << len(idx)
def cong(h, S):     # h S h^T over F_2; h as list of rows
    hS = [[sum(h[i][a] * S[a][j] for a in range(m)) & 1 for j in range(m)] for i in range(m)]
    return [[sum(hS[i][a] * h[j][a] for a in range(m)) & 1 for j in range(m)] for i in range(m)]
def ident(): return [[int(i == j) for j in range(m)] for i in range(m)]
gens = []
sw = ident(); sw[0], sw[1] = sw[1], sw[0]; gens.append(sw)
cyc = [[int(j == (i + 1) % m) for j in range(m)] for i in range(m)]; gens.append(cyc)
tv = ident(); tv[0][1] = 1; gens.append(tv)
G = [[enc(cong(g, dec(v))) for v in range(nS)] for g in gens]
orb = [-1] * nS; orbits = []
for v in range(nS):
    if orb[v] >= 0: continue
    k = len(orbits); orb[v] = k; st = [v]; mem = [v]
    while st:
        a = st.pop()
        for g in G:
            b = g[a]
            if orb[b] < 0: orb[b] = k; st.append(b); mem.append(b)
    orbits.append(mem)
def rank2(S):
    rows = [sum(S[i][j] << j for j in range(m)) for i in range(m)]; r = 0
    for c in range(m):
        p = next((i for i in range(r, m) if (rows[i] >> c) & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(m):
            if i != r and (rows[i] >> c) & 1: rows[i] ^= rows[r]
        r += 1
    return r
for k, mem in enumerate(orbits):
    S = dec(mem[0]); alt = all(S[i][i] == 0 for i in range(m))
    print("orbit %d: size %d, rank %d, %s" % (k, len(mem), rank2(S), "alternating" if alt else "non-alternating"))
I = enc(ident()); z = enc(dec(1))   # e1 e1^T
best = (0, None)
for mask in range(1, 1 << len(orbits)):
    U = set(v for k, mem in enumerate(orbits) if (mask >> k) & 1 for v in mem)
    fz = sum(1 for v in U if (v ^ z) not in U)
    fI = sum(1 for v in U if (v ^ I) not in U)
    if fz and fI / fz > best[0]:
        best = (fI / fz, mask, fI, fz)
print("best ratio f(I)/f(zz^T) over unions of orbits: %.6f (mask %s, f(I) = %s, f(move) = %s); plain needs %d" % (best[0], bin(best[1]), best[2], best[3], m))
