# fg_build.py -- every irreducible of Sp(6,2) as explicit real orthogonal matrices, cut out of the monomial modules
# Ind_B^G(lambda) (B = stabiliser of a flag inside L0, a Sylow 2-subgroup of P_0; 16 linear characters lambda), by the
# eigenspaces of the class sums of 2a (transvections) and 2b. Saves, per irreducible, the matrices of: the 6 generators of
# A = Sym(3,2), the 2 generators of the Levi GL(3,2), h, h2, k, and all 11 generators of G (for checks).
# Usage: fg_build.py out.npz [alpha]
import sys, json, itertools, time
import numpy as np
t0 = time.time()
m = 3; NV = 64
def Bf(u, v): return (bin((u & 7) & (v >> 3)).count('1') + bin((v & 7) & (u >> 3)).count('1')) & 1
def trans(w): return np.array([v ^ (w if Bf(v, w) else 0) for v in range(NV)], dtype=np.int64)
def comp(g, h): return g[h]                    # (g o h)(v) = g(h(v))
def inv(g): r = np.empty_like(g); r[g] = np.arange(NV); return r
def key(g): return g.astype(np.int8).tobytes()
def span(vs):
    S = {0}
    for v in vs: S |= {x ^ v for x in S}
    return frozenset(S)
X = lambda i: 1 << (i - 1); Z = lambda i: 1 << (2 + i)          # (e_i|0), (0|e_i)
EE = lambda i: X(i) | Z(i)
gens_w = [X(1), X(2), X(3), Z(1), Z(2), Z(3), X(1) | X(2), X(2) | X(3), Z(1) | Z(2), Z(2) | Z(3), EE(1)]
GENS = [trans(w) for w in gens_w]
A_GENS = [trans(w) for w in [Z(1), Z(2), Z(3), Z(1) | Z(2), Z(1) | Z(3), Z(2) | Z(3)]]
def levi(Amat):    # (x|z) -> (A x | A^{-T} z), Amat 3x3 over F_2 (list of rows)
    Am = np.array(Amat) % 2
    # inverse over F_2
    M = np.concatenate([Am, np.eye(3, dtype=int)], axis=1)
    for c in range(3):
        p = next(r for r in range(c, 3) if M[r, c]); M[[c, p]] = M[[p, c]]
        for r in range(3):
            if r != c and M[r, c]: M[r] ^= M[c]
    AinvT = M[:, 3:].T
    def ap(Mt, x): return sum(((int(Mt[i] @ np.array([(x >> j) & 1 for j in range(3)])) & 1) << i) for i in range(3))
    return np.array([ap(Am, v & 7) | (ap(AinvT, v >> 3) << 3) for v in range(NV)], dtype=np.int64)
L_GENS = [levi([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), levi([[0, 0, 1], [1, 0, 0], [0, 1, 0]])]
h1 = trans(EE(1)); h2 = comp(h1, trans(EE(2))); k3 = comp(h2, trans(EE(3)))
HS = [h1, h2, k3]
# ---- frames and flags
frames = set()
for vs in itertools.combinations(range(1, NV), 3):
    S = span(vs)
    if len(S) == 8 and all(Bf(a, b) == 0 for a in S for b in S): frames.add(S)
frames = sorted(frames, key=lambda S: sorted(S))
L0 = span([Z(1), Z(2), Z(3)])
flags = []
for L in frames:
    pts = sorted(L - {0})
    lines = set()
    for a, b in itertools.combinations(pts, 2): lines.add(span([a, b]))
    for ln in sorted(lines, key=sorted):
        for p in sorted(ln - {0}): flags.append((p, ln, L))
FIDX = {f: i for i, f in enumerate(flags)}; NF = len(flags); assert NF == 2835, NF
F0 = (Z(1), span([Z(1), Z(2)]), L0); i0 = FIDX[F0]
def img(g, f): p, ln, L = f; return (int(g[p]), frozenset(int(g[x]) for x in ln), frozenset(int(g[x]) for x in L))
def fperm(g): return np.array([FIDX[img(g, f)] for f in flags], dtype=np.int64)
GP = [fperm(g) for g in GENS]
# coset representatives gamma_F (gamma_F F0 = F) by BFS
gam = [None] * NF; gam[i0] = np.arange(NV); order = [i0]
for i in order:
    for g, p in zip(GENS, GP):
        j = p[i]
        if gam[j] is None: gam[j] = comp(g, gam[i]); order.append(j)
assert all(x is not None for x in gam)
GAMINV = [inv(x) for x in gam]
print(f"flags {NF}, coset reps done ({time.time()-t0:.1f}s)", flush=True)
# ---- B = Stab(F0): Schreier generators, closure, derived subgroup, abelianisation
sch = {}
for i in range(NF):
    for g, p in zip(GENS, GP):
        b = comp(GAMINV[p[i]], comp(g, gam[i])); sch[key(b)] = b
Bel = {key(np.arange(NV)): np.arange(NV)}; fr = list(Bel.values())
sg = list(sch.values())
while fr:
    nx = []
    for x in fr:
        for s in sg:
            y = comp(x, s); ky = key(y)
            if ky not in Bel: Bel[ky] = y; nx.append(y)
    fr = nx
assert len(Bel) == 512, len(Bel)
Bl = list(Bel.values())
comms = {}
for x in Bl:
    xi = inv(x)
    for y in Bl:
        c = comp(comp(xi, inv(y)), comp(x, y)); comms[key(c)] = c
D = {key(np.arange(NV)): np.arange(NV)}; fr = list(D.values()); cg = list(comms.values())
while fr:
    nx = []
    for x in fr:
        for s in cg:
            y = comp(x, s); ky = key(y)
            if ky not in D: D[ky] = y; nx.append(y)
    fr = nx
nq = 512 // len(D); print("B order 512, |B'| =", len(D), " |B/B'| =", nq, flush=True)
coset = {}; reps = []
for x in Bl:
    if key(x) in coset: continue
    c = len(reps); reps.append(x)
    for d in D.values(): coset[key(comp(x, d))] = c
# F_2 coordinates on B/B' (elementary abelian, checked)
def cmul(a, b): return coset[key(comp(reps[a], reps[b]))]
basis = []; spanq = {0: 0}            # coset -> coordinate bits
for c in range(nq):
    if c in spanq: continue
    bi = len(basis); basis.append(c)
    new = {}
    for s, bits in spanq.items(): new[cmul(s, c)] = bits | (1 << bi)
    spanq.update(new)
assert len(spanq) == nq and all(cmul(a, a) == 0 for a in range(nq))
for a in range(nq):
    for b in range(nq): assert spanq[cmul(a, b)] == spanq[a] ^ spanq[b]
kb = len(basis); COORD = {kk: spanq[c] for kk, c in coset.items()}
print("abelianisation F_2^%d" % kb, flush=True)
def mono(g, chi):
    """the monomial matrix of g on Ind_B^G(lambda_chi): returns (P, s) with g e_i = s_i e_{P_i}"""
    P = fperm(g)
    s = np.empty(NF)
    for i in range(NF):
        b = comp(GAMINV[P[i]], comp(g, gam[i]))
        s[i] = -1.0 if bin(COORD[key(b)] & chi).count('1') & 1 else 1.0
    return P, s
# class 2a = 63 transvections; class 2b = t_u t_v t_{u+v} over the 315 totally isotropic lines
TR = [trans(w) for w in range(1, NV)]
lines = set()
for u in range(1, NV):
    for v in range(u + 1, NV):
        if Bf(u, v) == 0: lines.add(frozenset({u, v, u ^ v}))
assert len(lines) == 315
C2B = []
for ln in sorted(lines, key=sorted):
    a, b, c = sorted(ln); C2B.append(comp(comp(trans(a), trans(b)), trans(c)))
ct = json.load(open('gap/ctbl.json'))
ia, ib = ct['classes'].index('2a'), ct['classes'].index('2b')
assert ct['centralizers'][ia] == 1451520 // 63 and ct['centralizers'][ib] == 1451520 // 315
TABLE = [(chi[0], chi[ia], chi[ib]) for chi in ct['irr']]
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6180339887
out = {}; found = {}
for chi in range(1 << kb):
    t1 = time.time()
    Zm = np.zeros((NF, NF))
    MA = [mono(g, chi) for g in TR]; MB = [mono(g, chi) for g in C2B]
    for (P, s) in MA: Zm[P, np.arange(NF)] += s
    Za = Zm.copy()
    for (P, s) in MB: Zm[P, np.arange(NF)] += alpha * s
    assert np.allclose(Zm, Zm.T)
    w, V = np.linalg.eigh(Zm)
    cl = []; st = 0
    for i in range(1, NF + 1):
        if i == NF or w[i] - w[i - 1] > 1e-6: cl.append((st, i)); st = i
    line = []
    for (a, b) in cl:
        Q = V[:, a:b]; dimc = b - a
        ra = np.trace(Q.T @ Za @ Q) / dimc          # = 63 chi(2a)/chi(1)
        rb = (w[a:b].mean() - ra) / alpha            # = 315 chi(2b)/chi(1)
        cand = [j for j, (d, va, vb) in enumerate(TABLE) if abs(63 * va / d - ra) < 1e-6 and abs(315 * vb / d - rb) < 1e-6 and dimc % d == 0]
        if len(cand) != 1: line.append(f"?{dimc}"); continue
        j = cand[0]; mult = dimc // TABLE[j][0]; line.append(f"{TABLE[j][0]}({TABLE[j][1]})x{mult}")
        if mult == 1 and j not in found:
            found[j] = chi
            def rest(g):
                P, s = mono(g, chi); MQ = np.zeros_like(Q); MQ[P] = s[:, None] * Q; return Q.T @ MQ
            mats = {'A': np.array([rest(g) for g in A_GENS]), 'L': np.array([rest(g) for g in L_GENS]),
                    'H': np.array([rest(g) for g in HS]), 'G': np.array([rest(g) for g in GENS])}
            for nm, arr in mats.items(): out[f"irr{j}_{nm}"] = arr
    print(f"lambda {chi:2d}: {' + '.join(line)}   ({time.time()-t1:.1f}s)", flush=True)
    if len(found) == 30: break
print("irreducibles found:", sorted(found), "missing:", [j for j in range(30) if j not in found], flush=True)
np.savez_compressed(sys.argv[1], **out, table=np.array(TABLE), found=np.array(sorted(found)))
print(f"total {time.time()-t0:.1f}s")
