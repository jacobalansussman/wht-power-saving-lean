# v-ex2-certnum (skeptic): independent exact check of ex2-certnum's Computations 12-14 FROM THE DUMPED EXACT DATA
# (v56/*.json: h, k, h2 in the type-adapted basis; dims of the P_0-types).  Own code, Fractions only.
#  - block ranks of h, k, h2; need_tau = rank k[tau, other], cap_tau = rank h[tau, other]
#  - for two-piece data: my own inputs-only construction (design differs from v56.py: right inverse on pivot columns),
#    then the exact containment test of Lemma 8b in the projected form, AND a full program evaluation with outputs solved.
# Usage: python3 vgame.py data.json [W M]
import sys, json
from fractions import Fraction as Fr
D = json.load(open(sys.argv[1])); dims = D['dims']; n = sum(dims)
h = [[Fr(x) for x in r] for r in D['h']]; k = [[Fr(x) for x in r] for r in D['k']]; h2 = [[Fr(x) for x in r] for r in D['h2']]
off = [sum(dims[:i]) for i in range(len(dims))]; idx = [list(range(off[i], off[i] + dims[i])) for i in range(len(dims))]
def sub(M, R, C): return [[M[r][c] for c in C] for r in R]
def rref(Mrows):
    M = [list(r) for r in Mrows]; rk = 0; piv = []
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; iv = 1 / M[rk][c]; M[rk] = [x * iv for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        piv.append(c); rk += 1
    return M[:rk], piv
def rank(M): return len(rref(M)[0]) if M and M[0] else 0
r = len(dims)
print('dims', dims, 'n', n)
for name, X in (('h', h), ('k', k), ('h2', h2)):
    print('block ranks of %s:' % name, [[rank(sub(X, idx[a], idx[b])) for b in range(r)] for a in range(r)])
oth = lambda a: [c for b in range(r) if b != a for c in idx[b]]
need = [rank(sub(k, idx[a], oth(a))) for a in range(r)]; need2 = [rank(sub(h2, idx[a], oth(a))) for a in range(r)]
cap = [rank(sub(h, idx[a], oth(a))) for a in range(r)]
print('need(k) =', need, ' need(h2) =', need2, ' cap =', cap)
print('bound per kernel array: max need/cap =', max(Fr(need[a], cap[a]) for a in range(r) if cap[a]) if all(cap[a] or not need[a] for a in range(r)) else 'inf')
# h invertible? and some identities
print('rank h =', rank(h), ' rank k =', rank(k))
if len(sys.argv) > 3:
    W, Mv = int(sys.argv[2]), int(sys.argv[3]); assert r == 2
    # sources as n x nW matrices (adapted basis both sides); inputs E_l
    def zeros(a, b): return [[Fr(0)] * b for _ in range(a)]
    srcs = []
    for l in range(W):
        E = zeros(n, n * W)
        for i in range(n): E[i][n * l + i] = Fr(1)
        srcs.append(E)
    # needed gain vectors for piece a: rows of k[a, b] placed in b-coords of block l (basis via rref)
    gains = []
    for a in range(2):
        b = 1 - a; Gb, _ = rref(sub(k, idx[a], idx[b]))
        vs = []
        for l in range(W):
            for g in Gb:
                v = [Fr(0)] * (n * W)
                for jj, c in enumerate(idx[b]): v[n * l + c] = g[jj]
                vs.append(v)
        gains.append(vs)
    # factor h[a, b] = C (Rrows), Rrows = independent rows; right inverse on pivot columns
    fac = []
    for a in range(2):
        b = 1 - a; H = sub(h, idx[a], idx[b]); Rr, pv = rref(H)    # rows of Rr span the row space of H
        fac.append((Rr, pv))
    ptr = [0, 0]; mixes = []
    for j in range(Mv):
        c = zeros(n, n * W)
        for a in range(2):
            b = 1 - a; Rr, pv = fac[a]; ra = len(Rr)
            V = gains[a][ptr[a]:ptr[a] + ra]; ptr[a] += len(V)
            # want Rr @ c_b = V (rows); Rr is in RREF with pivots pv (positions inside idx[b]): set c_b row pv[q] = V[q]
            for q, v in enumerate(V):
                row = idx[b][pv[q]]
                c[row] = list(v)
        # c must be a free combination of the inputs: c's a-rows supported on a-coords (any block) -- check
        for a in range(2):
            for i in idx[a]:
                assert all(c[i][n * l + cc] == 0 for l in range(W) for cc in range(n) if cc not in idx[a])
        # u = h c
        u = [[sum((h[i][q] * c[q][jj] for q in range(n) if h[i][q] != 0 and c[q][jj] != 0), Fr(0)) for jj in range(n * W)] for i in range(n)]
        srcs.append(u)
    print('gains transferred:', ptr, 'of', [len(g) for g in gains])
    # outputs: for each piece a, the a-rows of target (k in block t) must lie in span of the a-rows of all sources
    ok = True
    for a in range(2):
        stack = [srcs[s_][i] for s_ in range(len(srcs)) for i in idx[a]]
        base, _ = rref(stack); rb = len(base)
        for t in range(W):
            tg = []
            for i in idx[a]:
                v = [Fr(0)] * (n * W)
                for jj in range(n): v[n * t + jj] = k[i][jj]
                tg.append(v)
            r2 = rank(base + tg)
            print('piece %d (dim %d), output %d: rank sources %d, with target rows %d -> %s' % (a, dims[a], t, rb, r2, 'CONTAINED' if r2 == rb else 'MISSING'))
            ok &= r2 == rb
    print('VERDICT %d in %d (inputs-only, own design): %s' % (W, Mv, 'EXISTS (Lemma 8a/b)' if ok else 'not closed by this design'))
