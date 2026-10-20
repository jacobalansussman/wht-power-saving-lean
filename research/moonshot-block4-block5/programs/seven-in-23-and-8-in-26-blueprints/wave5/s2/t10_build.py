"""t10_build.py (session 2) -- build the ACTUAL O O F O F skeleton with a group-algebra middle content and decide the
full frame-sum relaxation with relax.Eng.decide (top content free; independent of the reduction in t9_galg.py).
Skeleton (Z7, 7 in 23): O1 = single-term start (cell Lam_1), O2 = single term from O1 (cell Lam_t), F1 = absorbed
shared move m . sum_h psi(h)^-1 beta_h tau_h(y) with y = sum_X y_X l_X t, letters l_X in Stab(L_Z) cap Stab(Lam_m),
i.e. F1's terms are  y_X psi(h)^-1 [array h, state c'_h^-1 o X o pt^-h o RF]  (X in Stab(Lam_t) cap Stab(Lam_F)),
F2 terminal (exempt, any character).  Y = sum_X y_X X is solved from the reduced problem; then decide()."""
import sys, os, random, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from relax import *
t0 = time.time()
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 21
psi = int(sys.argv[2]) if len(sys.argv) > 2 else 1
E = Eng('Z7', seed=seed); P = E.P; D = P.D; rng = random.Random(seed); nrng = np.random.default_rng(seed)
iZ, iK = E.iZ, E.iK
gensM = E.Dsig + E._letters[:30]
trans = {iZ: IDM}; fr = [iZ]
while fr:
    nf = []
    for i in fr:
        for g in gensM:
            s2 = mul(g, trans[i]); j = P.idx[span_mask([app(s2, 1 << (4 + k)) for k in range(4)])]
            if j not in trans: trans[j] = s2; nf.append(j)
    fr = nf
F = [j for j in P.nb[iK] if P.dZ[j] == 3][0]
nbF = P.nb[F]; ni = {Y: k for k, Y in enumerate(nbF)}
T = [Y for Y in nbF if P.dZ[Y] == 2 and D[Y, iK] == 2]
t1 = T[0]
def pm(a, b): return tuple(a[b[i]] for i in range(len(b)))
def pinv(a):
    r = [0] * len(a)
    for x, y in enumerate(a): r[y] = x
    return tuple(r)
def sp_pow(g, e):
    r = IDM
    for _ in range(e): r = mul(g, r)
    return r
def sp_order(g):
    k = 1; h = g
    while h != IDM: h = mul(g, h); k += 1
    return k
# group G (perm on the 30 slots -> a symplectic element), content cell t1
s0 = trans[t1]; s0i = inv(s0)
Fp = P.idx[img_mask(table(s0i), P.L[F])]
gens = []
for b in E._letters + [sigma(G.rand_engine(rng)) for _ in range(3000)]:
    if P.idx[img_mask(table(b), P.L[Fp])] == Fp:
        g = mul(s0, mul(b, s0i)); Tg = table(g)
        gens.append((tuple(ni[P.idx[img_mask(Tg, P.L[Y])]] for Y in nbF), g))
    if len(gens) >= 40: break
idp = tuple(range(30)); grp = {idp: IDM}; fr = [idp]
while fr:
    nf = []
    for g in fr:
        for (s, sg) in gens:
            h = pm(s, g)
            if h not in grp: grp[h] = mul(sg, grp[g]); nf.append(h)
    fr = nf
print('G order', len(grp), ' %.0fs' % (time.time() - t0), flush=True)
# check: the symplectic elements fix Lam_t and Lam_F and act on the slots as recorded
for g, gs in list(grp.items())[:50]:
    Tg = table(gs)
    assert P.idx[img_mask(Tg, P.L[t1])] == t1 and P.idx[img_mask(Tg, P.L[F])] == F
    assert tuple(ni[P.idx[img_mask(Tg, P.L[Y])]] for Y in nbF) == g
# twist pt: an order-7 element (exact order 7 as a symplectic matrix)
cands = []
for g, gs in grp.items():
    h = g; k = 1
    while h != idp and k < 8: h = pm(g, h); k += 1
    if k == 7: cands.append(g)
cands.sort(); ptp = cands[int(sys.argv[3]) if len(sys.argv) > 3 else 0]
gs = grp[ptp]; o = sp_order(gs); a2 = 0
while o % 2 == 0: o //= 2; a2 += 1
ept = (2 ** a2) * pow(2 ** a2, -1, 7); pt = sp_pow(gs, ept)
assert sp_order(pt) == 7 and tuple(ni[P.idx[img_mask(table(pt), P.L[Y])]] for Y in nbF) == ptp
# chain O1 -> O2 -> move m (as t6_pairbuild.py)
L1c = [Y for Y in P.nb[t1] if P.dZ[Y] == 1]
lam1 = L1c[int(sys.argv[4]) if len(sys.argv) > 4 else 0]
m1 = E.move_to(lam1); O1 = E.lmul(m1, E.inp(0)); (l, R1), = O1.keys()
lam_m2 = P.idx[img_mask(table(inv(R1)), P.L[t1])]
m2 = E.move_to(lam_m2); O2 = E.lmul(m2, O1); (l2, R2), = O2.keys()
assert E.cell(R2) == t1
lam_m = P.idx[img_mask(table(inv(R2)), P.L[F])]
m = E.move_to(lam_m)
RF = mul(R2, inv(m))
assert P.idx[span_mask([RF[4 + j] for j in range(4)])] == F
# basis of A from group elements (incremental Gram-Schmidt)
def pmat(g):
    M = np.zeros((30, 30)); M[list(g), list(range(30))] = 1; return M
keys = sorted(grp); random.Random(seed).shuffle(keys)
Q = np.zeros((0, 900)); basis = []
for g in keys:
    v = pmat(g).reshape(-1); r = v - Q.T @ (Q @ v)
    if np.linalg.norm(r) > 1e-8:
        Q = np.vstack([Q, r / np.linalg.norm(r)]); basis.append(g)
dA = len(basis)
print('dim A', dA, ' %.0fs' % (time.time() - t0), flush=True)
w7 = np.exp(2j * np.pi / 7)
def ovecs(g, a):
    gi = pinv(g); seen = set(); vs = []
    for N in range(30):
        if N in seen: continue
        orb = [N]; cur = N
        for _ in range(6): cur = gi[cur]; orb.append(cur)
        seen |= set(orb)
        if len(set(orb)) == 1: continue
        v = np.zeros(30, dtype=complex)
        for h, c in enumerate(orb): v[c] += w7 ** (a * h)
        vs.append(v)
    return np.array(vs).T
fixed = [s for s in range(30) if all(g[s] == s for (g, _) in gens)]
Fr_slots = sorted(set([ni[c] for c in nbF if D[c, lam1] == 1]) | set(fixed))
k0 = ni[iK]; e = np.zeros(30, dtype=complex); e[k0] = 1
Xi = [int(c) for c in os.environ.get('XI', '1,2,3,4,5,6').split(',')]   # effective characters a = chi - psi (default all six nontrivial: E = {psi})
Bm = [pmat(g) for g in basis]
nF = len(Fr_slots)
Mbig = np.zeros((30 * len(Xi), dA + nF * len(Xi)), dtype=complex); bb = np.zeros(30 * len(Xi), dtype=complex)
for ia, a in enumerate(Xi):
    S = ovecs(ptp, a) if a else np.concatenate([ovecs(ptp, 0), np.eye(30)[:, fixed]], axis=1); wa = S @ (nrng.standard_normal(S.shape[1]) + 1j * nrng.standard_normal(S.shape[1]))
    for k, B in enumerate(Bm): Mbig[30 * ia:30 * ia + 30, k] = B @ wa
    for j, s in enumerate(Fr_slots): Mbig[30 * ia + s, dA + nF * ia + j] = 1
    bb[30 * ia:30 * ia + 30] = e
x, *_ = np.linalg.lstsq(Mbig, bb, rcond=None)
print('reduced problem residual %.2e' % np.linalg.norm(Mbig @ x - bb), flush=True)
y = x[:dA]
# the actual shared move F1
def F1_of(coefs, elems):
    u = {}
    for c, g in zip(coefs, elems):
        if abs(c) < 1e-13: continue
        X = grp[g]; cur = IDM
        for h in range(7):
            Rh = mul(E.cpi[h], mul(X, mul(cur, RF)))
            k = (h, Rh); u[k] = u.get(k, 0) + c * np.conj(E.chars[psi][h])
            cur = mul(inv(pt), cur)
    return {k: v for k, v in u.items() if abs(v) > 1e-12}
F1 = F1_of(y, basis)
print('F1 terms', len(F1), ' max distance', E.maxdist(F1), ' %.0fs' % (time.time() - t0), flush=True)
ok, info = E.decide([E.inp(0), O1, O2], [F1], E=[psi])
print('DESIGNED (Y from the reduced problem): passes =', ok, ' mu_ok', info.get('mu_ok'), ' null_dim', info.get('null_dim'),
      ' sv_gap', info.get('sv_gap'), flush=True)
# controls: the single-term F1 (Y = identity), and a random Y with the same support
okI, infoI = E.decide([E.inp(0), O1, O2], [F1_of([1.0], [idp])], E=[psi])
print('CONTROL single term (Y = 1): passes =', okI, ' mu_ok', infoI.get('mu_ok'), flush=True)
yr = nrng.standard_normal(dA) + 1j * nrng.standard_normal(dA)
okR, infoR = E.decide([E.inp(0), O1, O2], [F1_of(yr, basis)], E=[psi])
print('CONTROL random Y (same support): passes =', okR, ' mu_ok', infoR.get('mu_ok'), flush=True)
np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'y_seed{seed}_psi{psi}.npy'), np.array(y))
print('done %.0fs' % (time.time() - t0))
