"""t14_z2cube_allsubs.py (session 2) -- EXHAUSTIVE exact test of the group-algebra family for 8 in 26 (Z2^3), per content
cell: every elementary abelian subgroup Q of order 8 of G (image of Stab(Lam_t) cap Stab(Lam_F) on the 30 slots) whose orbit
through the kernel slot is regular.  For a twist q with image Q, the top's F1 part in character a is Y W_a with Y in A
(span of G's permutation matrices) and W_a in S_a (sign-weighted Q-orbit vectors), so  e_K in A.S_a + F  is NECESSARY for
every a outside E (|E| <= 2).  Counts the characters a = 1..7 failing it, with F = one O1 cell (best of 3) and with the
generous F = every slot at input distance <= 3 (the most any O1 / O2 contents give on nbrs(Lam_F)).  A subgroup with >= 3
failing characters cannot give a unit in this family.  (Subgroups are taken at the permutation level: a superset of the
exactly liftable ones, so an exclusion here covers them.)"""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp8 import *
from collections import Counter
t0 = time.time()
ct = int(sys.argv[1])
P = Polar(); D = P.D; iZ = P.idx[LZ]
iK = P.idx[img_mask(table(inv(sigma(G.K))), LZ)]
rng = random.Random(3)
gensM = [sigma(G.D[z]) for z in range(1, 16)] + [sigma(G.rand_engine(rng)) for _ in range(30)]
trans = {iZ: IDM}; fr = [iZ]
while fr:
    nf = []
    for i in fr:
        for g in gensM:
            s2 = mul(g, trans[i]); j = P.idx[span_mask([app(s2, 1 << (4 + k)) for k in range(4)])]
            if j not in trans: trans[j] = s2; nf.append(j)
    fr = nf
letters = [sigma(G.rand_engine(rng)) for _ in range(4000)]
F = [j for j in P.nb[iK] if P.dZ[j] == 3][0]
nbF = P.nb[F]; ni = {Y: k for k, Y in enumerate(nbF)}
T = [Y for Y in nbF if P.dZ[Y] == 2 and D[Y, iK] == 2]
t1 = T[ct]
def pm(a, b): return tuple(a[b[i]] for i in range(len(b)))
s0 = trans[t1]; s0i = inv(s0)
Fp = P.idx[img_mask(table(s0i), P.L[F])]
gens = []
for b in letters:
    if P.idx[img_mask(table(b), P.L[Fp])] == Fp:
        g = mul(s0, mul(b, s0i)); Tg = table(g)
        gens.append(tuple(ni[P.idx[img_mask(Tg, P.L[Y])]] for Y in nbF))
    if len(gens) >= 60: break
idp = tuple(range(30)); grp = {idp}; fr = [idp]
while fr:
    nf = []
    for g in fr:
        for s in gens:
            h = pm(s, g)
            if h not in grp: grp.add(h); nf.append(h)
    fr = nf
k0 = ni[iK]
def pmat(g):
    M = np.zeros((30, 30)); M[list(g), list(range(30))] = 1; return M
keys = sorted(grp)
Q = np.zeros((0, 900)); basis = []
for g in keys:
    v = pmat(g).reshape(-1); r = v - Q.T @ (Q @ v)
    if np.linalg.norm(r) > 1e-8: Q = np.vstack([Q, r / np.linalg.norm(r)]); basis.append(g)
Bm = np.array([pmat(g) for g in basis]); dA = len(basis)
fixed = [s for s in range(30) if all(g[s] == s for g in gens)]
L1c = [Y for Y in P.nb[t1] if P.dZ[Y] == 1]
Fmain = [sorted(set([ni[c] for c in nbF if D[c, lam1] == 1]) | set(fixed)) for lam1 in L1c]
Fgen = [ni[c] for c in nbF if P.dZ[c] <= 3]
print('content cell', t1, ': |G|', len(grp), ' dim A', dA, ' %.0fs' % (time.time() - t0), flush=True)
invs = [g for g in grp if pm(g, g) == idp and g != idp]
comm = {}
for a in invs:
    for b in invs:
        if a < b and pm(a, b) == pm(b, a) and pm(a, b) != idp:
            comm.setdefault(a, []).append(b)
subs = set()
for a, bl in comm.items():
    for b in bl:
        ab = pm(a, b)
        for c in comm.get(b, []):
            if c in (a, ab) or pm(a, c) != pm(c, a): continue
            Esub = frozenset([idp, a, b, ab, c, pm(a, c), pm(b, c), pm(ab, c)])
            if len(Esub) == 8: subs.add(Esub)
regs = [S for S in subs if len({g[k0] for g in S}) == 8]
print('involutions', len(invs), '; elementary abelian subgroups of order 8:', len(subs), '; with a regular orbit through the kernel slot:', len(regs), ' %.0fs' % (time.time() - t0), flush=True)
e = np.zeros(30); e[k0] = 1
def inside(Mx):
    x, *_ = np.linalg.lstsq(Mx, e, rcond=None); return np.linalg.norm(Mx @ x - e) < 1e-8
cnt = Counter(); best = None
for n_, S in enumerate(regs):
    els = sorted(S)
    # a basis of the subgroup: pick 3 generators, label h by 3 bits
    g1 = els[1]; g2 = next(g for g in els if g not in (idp, g1)); g3 = next(g for g in els if g not in (idp, g1, g2, pm(g1, g2)))
    def qp(h):
        g = idp
        for bit, a in ((4, g1), (2, g2), (1, g3)):
            if h & bit: g = pm(a, g)
        return g
    qs = [qp(h) for h in range(8)]
    fm = 0; fg = 0
    for a in range(1, 8):
        sn = set(); vs = []
        for s_ in range(30):
            if s_ in sn: continue
            v = np.zeros(30)
            for h in range(8):
                c = qs[h][s_]; v[c] += (-1) ** bin(a & h).count('1'); sn.add(c)
            if np.linalg.norm(v) > 1e-9: vs.append(v)
        AS = np.einsum('kij,vj->ikv', Bm, np.array(vs)).reshape(30, -1)
        if not any(inside(np.concatenate([AS, np.eye(30)[:, Fs]], axis=1)) for Fs in Fmain): fm += 1
        if not inside(np.concatenate([AS, np.eye(30)[:, Fgen]], axis=1)): fg += 1
    cnt[(fm, fg)] += 1
    if best is None or fg < best[0]: best = (fg, fm, n_)
    if n_ % 200 == 0: print('  ', n_, 'of', len(regs), dict(cnt), ' %.0fs' % (time.time() - t0), flush=True)
print('content cell', t1, ': (#failing characters with one O1 cell, with every distance<=3 slot free) ->', dict(sorted(cnt.items())))
print('fewest failing characters (generous free span):', best[0] if best else None, '; EXCLUDED (>= 3 failing for every subgroup):', (best[0] >= 3) if best else None)
print('done %.0fs' % (time.time() - t0), flush=True)
