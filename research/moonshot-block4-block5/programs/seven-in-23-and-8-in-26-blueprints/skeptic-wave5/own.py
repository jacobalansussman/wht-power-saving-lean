"""own.py (v-w5-blueprint, skeptic): independent rebuild of the w5-blueprint reduced problem from F_2^8 symplectic data only.
No Clifford matrices, no sp8.py.  Lagrangians from transvection BFS; G = Stab(t) cap Stab(F) by Schreier generators;
A = span of G's 30x30 permutation matrices; Z7 test e_K in Y S_a + Fr for all a = 1..6 with ONE Y in A; controls;
Z2^3 test e_K in A.S_a + Fgen over ALL elementary abelian subgroups of order 8 (regular kernel orbit or not)."""
import sys, random, itertools, time
import numpy as np
t0 = time.time()
POP = np.array([bin(i).count('1') for i in range(256)])
def om(u, v): return (POP[(u & 15) & (v >> 4)] + POP[(u >> 4) & (v & 15)]) & 1
V = np.arange(256)
OM = np.array([[om(u, v) for v in range(256)] for u in range(256)])
def transv(v): return np.where(OM[:, v] == 1, V ^ v, V).astype(np.int64)   # T_v(x) = x + om(x,v) v
TR = [transv(v) for v in range(1, 256)]
def mask(S): 
    m = 0
    for x in S: m |= 1 << int(x)
    return m
def elems(mk): return [i for i in range(256) if (mk >> i) & 1]
def img(g, mk): return mask(g[elems(mk)])
LZ = mask([v for v in range(256) if v & 15 == 0])
LX = mask([v for v in range(256) if v >> 4 == 0])
# BFS all Lagrangians, with transversal elements tr[L] mapping LZ -> L
idx = {LZ: 0}; Ls = [LZ]; tr = {LZ: V.copy()}; fr = [LZ]
while fr:
    nf = []
    for x in fr:
        for T in TR:
            y = img(T, x)
            if y not in idx:
                idx[y] = len(Ls); Ls.append(y); tr[y] = T[tr[x]]; nf.append(y)
    fr = nf
assert len(Ls) == 2295, len(Ls)
def dim(mk): return bin(mk).count('1').bit_length() - 1
def dist(a, b): return 4 - dim(a & b)
# check isotropy of all
for x in Ls[:50]:
    el = elems(x); assert all(OM[u, v] == 0 for u in el for v in el)
print('Lagrangians', len(Ls), '%.0fs' % (time.time() - t0), flush=True)
Lset = Ls
def nbrs(x): return [y for y in Lset if dist(x, y) == 1]
K = LX; assert dist(K, LZ) == 4
nbK = nbrs(K)
Fc = [y for y in nbK if dist(y, LZ) == 3]
iF = int(sys.argv[1]) if len(sys.argv) > 1 else 0
F = Fc[iF]
nbF = nbrs(F); assert len(nbF) == 30
sl = {y: k for k, y in enumerate(nbF)}
k0 = sl[K]
Tc = [y for y in nbF if dist(y, LZ) == 2 and dist(y, K) == 2]
print('neighbours of K at d3:', len(Fc), '; content cells t:', len(Tc), ' %.0fs' % (time.time() - t0), flush=True)
def compose(a, b): return a[b]
def inverse(a):
    r = np.empty(256, dtype=np.int64); r[a] = V; return r
rng = random.Random(11)
def rand_sp(n=40):
    g = V.copy()
    for _ in range(n): g = TR[rng.randrange(255)][g]
    return g
# transversal of Sp on Lagrangians starting at F: trF[L] maps F -> L
trF = {}
gF = tr[F]   # LZ -> F
for y in Ls:
    trF[y] = compose(tr[y], inverse(gF))     # F -> LZ -> y
def stabF():
    g = rand_sp(); y = img(g, F); return compose(inverse(trF[y]), g)
def slotperm(g): return tuple(sl[img(g, y)] for y in nbF)
def build_G(t):
    # Schreier: random elements of Stab(F), then transversal on nbrs(F) orbit of t
    SF = [stabF() for _ in range(25)]
    for s in SF: assert img(s, F) == F
    # transversal of <SF> on nbF orbit of t
    trt = {t: V.copy()}; frn = [t]
    while frn:
        nf = []
        for x in frn:
            for s in SF:
                y = img(s, x)
                if y not in trt: trt[y] = compose(s, trt[x]); nf.append(y)
        frn = nf
    gens = []
    for _ in range(40):
        s = stabF(); y = img(s, t)
        if y not in trt: continue
        u = compose(inverse(trt[y]), s); assert img(u, t) == t and img(u, F) == F
        gens.append(u)
    pg = [slotperm(u) for u in gens]
    idp = tuple(range(30)); grp = {idp}; frn = [idp]
    while frn:
        nf = []
        for g in frn:
            for s in pg:
                h = tuple(s[g[i]] for i in range(30))
                if h not in grp: grp.add(h); nf.append(h)
        frn = nf
    return gens, sorted(grp)
def pmat(g):
    M = np.zeros((30, 30)); M[list(g), list(range(30))] = 1; return M
def algebra(grp):
    Q = np.zeros((0, 900)); basis = []
    for g in grp:
        v = pmat(g).reshape(-1); r = v - Q.T @ (Q @ v)
        if np.linalg.norm(r) > 1e-8: Q = np.vstack([Q, r / np.linalg.norm(r)]); basis.append(g)
    return np.array([pmat(g) for g in basis])

# ------------------------------------------------------------------ main
mode = sys.argv[2] if len(sys.argv) > 2 else 'z7'
it = int(sys.argv[3]) if len(sys.argv) > 3 else 0
t = Tc[it]
gens, grp = build_G(t)
Bm = algebra(grp); dA = len(Bm)
fixed = [s for s in range(30) if all(g[s] == s for g in grp)]
print(f'F choice {iF}, content cell #{it}: |G| = {len(grp)}, dim A = {dA}, G-fixed slots {fixed} (t slot {sl[t]}), input distances of fixed',
      [dist(nbF[s], LZ) for s in fixed], ' %.0fs' % (time.time() - t0), flush=True)
# G-orbits on slots
orbs = []; seen = set()
for s in range(30):
    if s in seen: continue
    o = sorted({g[s] for g in grp}); seen |= set(o); orbs.append(o)
print('G-orbit sizes on the 30 slots:', [len(o) for o in orbs], '; kernel slot', k0, 'in orbit of size', len([o for o in orbs if k0 in o][0]), flush=True)
L1c = [y for y in nbrs(t) if dist(y, LZ) == 1]
e = np.zeros(30, dtype=complex); e[k0] = 1
def pm(a, b): return tuple(a[b[i]] for i in range(30))
def pinv(a):
    r = [0] * 30
    for x, y in enumerate(a): r[y] = x
    return tuple(r)
def inspan(M, v, tol=1e-9):
    x, *_ = np.linalg.lstsq(M, v, rcond=None); return float(np.linalg.norm(M @ x - v))
import os
def free_slots(lam1):
    if os.environ.get('NOFREE'): return []
    return sorted({s for s in range(30) if dist(nbF[s], lam1) == 1} | {sl[t]} | {s for s in range(30) if dist(nbF[s], t) == 1})
if mode == 'z7':
    w = np.exp(2j * np.pi / 7)
    sev = []
    for g in grp:
        h = g; k = 1
        while h != tuple(range(30)) and k < 8: h = pm(g, h); k += 1
        if k == 7: sev.append(g)
    cls = {}; reps = []
    for g in sev:
        if g in cls: continue
        reps.append(g)
        for x in grp: cls.setdefault(pm(x, pm(g, pinv(x))), len(reps) - 1)
    print('order-7 elements', len(sev), '; G-classes', len(reps), [sum(1 for v in cls.values() if v == c) for c in range(len(reps))], flush=True)
    def Sa(pt, a):
        seen = set(); vs = []
        for N in range(30):
            if N in seen: continue
            orb = [N]; cur = N
            for _ in range(6): cur = pt[cur]; orb.append(cur)
            seen |= set(orb)
            if len(set(orb)) == 1: continue
            v = np.zeros(30, dtype=complex)
            for h, c in enumerate(orb): v[c] += w ** (a * h)
            vs.append(v)
        return np.array(vs).T
    nrng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    for lam1 in L1c:
        Fr = free_slots(lam1)
        for ci, pt in enumerate(reps):
            for chars in ([1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]):
                nF = len(Fr); na = len(chars)
                M = np.zeros((30 * na, dA + nF * na), dtype=complex); b = np.zeros(30 * na, dtype=complex)
                Ws = []
                for ia, a in enumerate(chars):
                    S = Sa(pt, a); wa = S @ (nrng.standard_normal(S.shape[1]) + 1j * nrng.standard_normal(S.shape[1])); Ws.append(wa)
                    M[30 * ia:30 * ia + 30, :dA] = np.einsum('kij,j->ik', Bm, wa)
                    for j, s in enumerate(Fr): M[30 * ia + s, dA + nF * ia + j] = 1
                    b[30 * ia:30 * ia + 30] = e
                x, *_ = np.linalg.lstsq(M, b, rcond=None); res = float(np.linalg.norm(M @ x - b))
                sv = np.linalg.svd(M, compute_uv=False)
                line = f'  lam1 {idx[lam1]} class {ci} chars {chars}: eqs {M.shape[0]} unknowns {M.shape[1]} rank {(sv > 1e-9 * sv[0]).sum()} residual {res:.2e}'
                if chars[0] == 1:
                    # controls: Y = I (single term) and random Y in A, W_a free: e_K in Y S_a + Fr ?
                    for name, Y in (('Y=I', np.eye(30)), ('Y random in A', np.einsum('k,kij->ij', nrng.standard_normal(dA), Bm))):
                        fails = 0
                        for a in chars:
                            S = Sa(pt, a); Mx = np.concatenate([Y @ S, np.eye(30)[:, Fr]], axis=1)
                            if inspan(Mx, e) > 1e-8: fails += 1
                        line += f' | control {name}: fails {fails}/6'
                    # witness check: Y from x, verify per character with Y fixed and W_a free
                    Ysol = np.einsum('k,kij->ij', x[:dA], Bm); worst = 0
                    for ia, a in enumerate(chars):
                        S = Sa(pt, a); Mx = np.concatenate([Ysol @ S, np.eye(30)[:, Fr]], axis=1); worst = max(worst, inspan(Mx, e))
                    line += f' | solved Y re-check worst {worst:.1e}'
                print(line, flush=True)
    print('done %.0fs' % (time.time() - t0))
if mode == 'z2':
    from collections import Counter
    idp = tuple(range(30))
    invs = [g for g in grp if pm(g, g) == idp and g != idp]
    invset = set(invs)
    comm = {a: [b for b in invs if b != a and pm(a, b) == pm(b, a)] for a in invs}
    subs = set()
    for a in invs:
        for b in comm[a]:
            if b < a: continue
            ab = pm(a, b)
            for c in comm[b]:
                if c in (a, ab) or pm(a, c) != pm(c, a): continue
                S = frozenset([idp, a, b, ab, c, pm(a, c), pm(b, c), pm(ab, c)])
                if len(S) == 8: subs.add(S)
    print('involutions', len(invs), '; elementary abelian subgroups of order 8:', len(subs), ' %.0fs' % (time.time() - t0), flush=True)
    Gk0 = [o for o in orbs if k0 in o][0]
    D4 = [s for s in range(30) if dist(nbF[s], LZ) == 4]
    lam1 = L1c[0]; Fr1 = free_slots(lam1)
    NF1 = [s for s in range(30) if s not in Fr1]
    print('distance-4 slots', len(D4), '; non-free slots with one O1 cell', len(NF1), flush=True)
    C = np.einsum('kis->sik', Bm)      # C[s] = 30 x dA, columns B_k e_s
    def fails_for(S):
        els = sorted(S)
        g1 = els[1]; g2 = next(g for g in els if g not in (idp, g1)); g3 = next(g for g in els if g not in (idp, g1, g2, pm(g1, g2)))
        qs = []
        for h in range(8):
            g = idp
            for bit, x in ((1, g1), (2, g2), (4, g3)):
                if h & bit: g = pm(x, g)
            qs.append(g)
        out = []
        for rows in (D4, NF1):
            nfail = 0; failing = []
            for a in range(1, 8):
                seen = set(); vs = []
                for s in range(30):
                    if s in seen: continue
                    v = np.zeros(30)
                    for h in range(8):
                        c = qs[h][s]; v[c] += (-1) ** (bin(a & h).count('1')); seen.add(c)
                    if np.abs(v).sum() > 1e-9: vs.append(v)
                if not vs: nfail += 1; failing.append(a); continue
                AS = np.concatenate([np.tensordot(v, C, axes=(0, 0)) for v in vs], axis=1)[rows]
                tgt = np.zeros(len(rows)); tgt[rows.index(k0)] = 1
                U, sv, _ = np.linalg.svd(AS, full_matrices=False); U = U[:, sv > 1e-9 * max(sv[0], 1e-30)]
                r = tgt - U @ (U.T @ tgt)
                if np.linalg.norm(r) > 1e-8: nfail += 1; failing.append(a)
            out.append(nfail)
        regk = len({g[k0] for g in S}) == 8
        regG = any(len({g[s] for g in S}) == 8 for s in Gk0)
        return out, regk, regG
    cnt = Counter(); minfail = {}
    for n_, S in enumerate(sorted(subs, key=lambda s: sorted(s))):
        (fg, f1), regk, regG = fails_for(S)
        cnt[(regk, regG, fg, f1)] += 1
        key = (regk, regG); minfail[key] = min(minfail.get(key, 99), fg)
        if n_ % 2000 == 0: print('  ', n_, ' %.0fs' % (time.time() - t0), flush=True)
    print('counts (regular kernel orbit, some regular orbit in G.k0, #failing chars with all d<=3 free, #failing with one O1 cell):')
    for k in sorted(cnt): print('   ', k, cnt[k])
    print('min failing (generous free span) by (regk, regG):', minfail)
    print('done %.0fs' % (time.time() - t0))
if mode == 'cal':
    # bilinear reduced problem: exists Y in A, w_a in S_a (a in chars), f_a in Fr with Y w_a + f_a = e_K.  ALS from random starts.
    w7 = np.exp(2j * np.pi / 7)
    sev = [g for g in grp if g != tuple(range(30)) and pm(g, pm(g, pm(g, pm(g, pm(g, pm(g, g)))))) == tuple(range(30))]
    pt = sev[0]
    def Sa(pt, a):
        seen = set(); vs = []
        for N in range(30):
            if N in seen: continue
            orb = [N]; cur = N
            for _ in range(6): cur = pt[cur]; orb.append(cur)
            seen |= set(orb)
            if len(set(orb)) == 1: continue
            v = np.zeros(30, dtype=complex)
            for h, c in enumerate(orb): v[c] += w7 ** (a * h)
            vs.append(v)
        return np.array(vs).T
    nrng = np.random.default_rng(5)
    for chars in ([1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]):
        for lam1 in L1c[:1]:
            Fr = free_slots(lam1); IF = np.eye(30)[:, Fr]
            Ss = [Sa(pt, a) for a in chars]
            best = 9
            for start in range(8):
                ws = [S @ (nrng.standard_normal(S.shape[1]) + 1j * nrng.standard_normal(S.shape[1])) for S in Ss]
                hist = []
                for itn in range(150):
                    # Y step
                    na = len(chars); nF = len(Fr)
                    M = np.zeros((30 * na, dA + nF * na), dtype=complex); b = np.tile(e, na)
                    for ia, wa in enumerate(ws):
                        M[30 * ia:30 * ia + 30, :dA] = np.einsum('kij,j->ik', Bm, wa); M[30 * ia:30 * ia + 30, dA + nF * ia: dA + nF * (ia + 1)] = IF
                    x, *_ = np.linalg.lstsq(M, b, rcond=None); res = np.linalg.norm(M @ x - b)
                    Y = np.einsum('k,kij->ij', x[:dA], Bm)
                    # w step (normalise Y to avoid drift)
                    nY = np.linalg.norm(Y); Y = Y / nY
                    ws = []; r2 = 0
                    for S in Ss:
                        Mx = np.concatenate([Y @ S, IF], axis=1); c, *_ = np.linalg.lstsq(Mx, e, rcond=None)
                        r2 += np.linalg.norm(Mx @ c - e) ** 2; ws.append(S @ c[:S.shape[1]])
                    hist.append(res)
                best = min(best, min(hist))
                print(f'  chars {chars} lam1 {idx[lam1]} start {start}: ALS residual first {hist[0]:.3e} -> last {hist[-1]:.3e} (min {min(hist):.3e})', flush=True)
            print(f' chars {chars}: best ALS residual {best:.3e}', flush=True)
    print('done %.0fs' % (time.time() - t0))
if mode == 'a0':
    w7 = np.exp(2j * np.pi / 7)
    sev = [g for g in grp if g != tuple(range(30)) and pm(g, pm(g, pm(g, pm(g, pm(g, pm(g, g)))))) == tuple(range(30))]
    for pt in sev[:3]:
        for lam1 in L1c:
            Fr = free_slots(lam1)
            for a in range(7):
                seen = set(); vs = []
                for N in range(30):
                    if N in seen: continue
                    orb = [N]; cur = N
                    for _ in range(6): cur = pt[cur]; orb.append(cur)
                    seen |= set(orb)
                    v = np.zeros(30, dtype=complex)
                    for h, c in enumerate(orb): v[c] += w7 ** (a * h)
                    if np.abs(v).sum() > 1e-9: vs.append(v)
                AS = np.concatenate([np.einsum('kij,j->ik', Bm, v) for v in vs] + [np.eye(30)[:, Fr]], axis=1)
                print(f'  pt {sev.index(pt)} lam1 {idx[lam1]} a {a}: dist(e_K, A.S_a + Fr) = {inspan(AS, e):.4f}', flush=True)
