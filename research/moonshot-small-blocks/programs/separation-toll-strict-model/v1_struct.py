# v1_struct.py -- v-ms-round2: recompute the key numbers of Theorem S8 with own code.
#  (a) dims N, E_R, E_K; failure ranks of D1, D2 on E_R; dim K_D1 cap K_D2
#  (b) the FULL Levi group of R = D1 D2 in the widened game (all monomials with Z4 phases h, R h R^-1 monomial),
#      enumerated directly (perms must preserve the cosets of <1,2>; all phases), not as diag x perm
#  (c) the orbit of K_D1 under the full Levi group; pairwise / triple / 4-way intersection dims (the grid)
import sys, time, pickle, itertools, numpy as np
from vlib import *

p = int(sys.argv[1]) if len(sys.argv) > 1 else P_DEFAULT
F = Fp(p); t0 = time.time()
print(f"p = {p}, i = {F.i}")
D1, D2, D4 = dirz(F, 1), dirz(F, 2), dirz(F, 4)
# homomorphism sanity
g1 = mmod(F, D1, mono(F, [3, 1, 0, 2, 7, 5, 6, 4], [1, 0, 2, 3, 0, 1, 1, 2]))
g2 = mmod(F, D4, D2)
assert np.array_equal(rho(F, mmod(F, g1, g2)), F.mm(rho(F, g1), rho(F, g2))), "rho not multiplicative"
print("rho multiplicative on a test pair: OK")
R = mmod(F, D1, D2); K = mmod(F, D1, D2, D4)
rR = rho(F, R); rK = rho(F, K)
print("dim N =", len(NIDX))
ER = F.null(rR[np.ix_(OIDX, NIDX)])            # N coords, 280 x dimE
EK = F.null(rK[np.ix_(OIDX, NIDX)])
print("dim E_R =", ER.shape[1], " dim E_K =", EK.shape[1], " demands", 280 - ER.shape[1], 280 - EK.shape[1])
rD1 = rho(F, D1); rD2 = rho(F, D2); rD4 = rho(F, D4)
print("cap: rank pi_O rho(D1)|_N =", F.rank(rD1[np.ix_(OIDX, NIDX)]))
A1 = F.mm(rD1[np.ix_(OIDX, NIDX)], ER); A2 = F.mm(rD2[np.ix_(OIDX, NIDX)], ER); A4 = F.mm(rD4[np.ix_(OIDX, NIDX)], ER)
print("failure ranks on E_R: D1", F.rank(A1), " D2", F.rank(A2), " D4", F.rank(A4))
K1 = F.mm(ER, F.null(A1)); K2 = F.mm(ER, F.null(A2))   # N coords
print("dim K_D1", K1.shape[1], " dim K_D2", K2.shape[1], " dim K_D1 cap K_D2 =",
      K1.shape[1] + K2.shape[1] - F.rank(np.hstack([K1, K2])))

# ---- (b) full Levi group (complex floats; entries Gaussian dyadic so tolerance is safe)
def cdir(z):
    M = np.zeros((8, 8), complex)
    for x in range(8): M[x, x] += (1 + 1j) / 2; M[x ^ z, x] += (1 - 1j) / 2
    return M
Rc = cdir(1) @ cdir(2); RcH = Rc.conj().T
blocks = [[0, 1, 2, 3], [4, 5, 6, 7]]
perms = []
for sw in [0, 1]:
    for a in itertools.permutations(range(4)):
        for b in itertools.permutations(range(4)):
            pm = [0] * 8
            for x in range(4): pm[x] = blocks[sw][a[x]]; pm[4 + x] = blocks[1 - sw][b[x]]
            perms.append(pm)
# sanity: a perm NOT preserving the blocks never works (sample)
phs = np.array([ph for ph in itertools.product(range(4), repeat=8) if ph[0] == 0])     # 16384 x 8
phase = (1j) ** phs                                                                     # 16384 x 8
def levi_phases(pm):
    Pm = np.zeros((8, 8), complex)
    for x in range(8): Pm[pm[x], x] = 1
    A = Rc @ Pm                                     # (R P)[a,k]
    # (R P diag(phase) R^H)[n,a,b] = sum_k A[a,k] phase[n,k] RcH[k,b]
    Mt = np.einsum('ak,nk,kb->nab', A, phase, RcH, optimize=True)
    nzc = (np.abs(Mt) > 1e-9).sum(axis=2)           # nonzeros per row
    ok = (nzc == 1).all(axis=1)
    return np.nonzero(ok)[0]
bad = 0
rng = np.random.default_rng(5)
for _ in range(40):
    pm = list(rng.permutation(8))
    if any(sorted(pm[x] for x in range(4)) not in ([0, 1, 2, 3], [4, 5, 6, 7]) for _ in [0]) and \
       sorted(pm[x] for x in range(4)) not in ([0, 1, 2, 3], [4, 5, 6, 7]):
        if len(levi_phases(pm)): bad += 1
print("non-block-preserving perms with a Levi phase (sample of up to 40):", bad)
levi = []
for pm in perms:
    for n in levi_phases(pm): levi.append((tuple(pm), tuple(int(v) for v in phs[n])))
print(f"FULL Levi group mod global phase: {len(levi)} elements  ({time.time()-t0:.0f}s)")
nd = sum(1 for pm, ph in levi if pm == tuple(range(8)))
npm = sum(1 for pm, ph in levi if all(v == 0 for v in ph))
print("  diagonal ones:", nd, "  pure permutations:", npm, "  distinct perms occurring:", len(set(pm for pm, ph in levi)))
pickle.dump(levi, open('v_levi.pkl', 'wb'))

# ---- (c) orbit of K_D1 under the full Levi group (rho(h) K_D1 for every h)
def ann(Kb):                       # rows whose common kernel (in N coords) is span(Kb)
    return F.null(Kb.T).T
posN = {int(j): k for k, j in enumerate(NIDX)}
def rhoN(pm, ph, X):               # rho(h) X for X in N coords (280 x c)
    tgt, ex = rho_mono_index(pm, ph)
    full = np.zeros((560, X.shape[1])); full[NIDX] = X
    out = apply_mono(F, tgt, ex, full)
    assert not out[OIDX].any()
    return out[NIDX]
orbit = [K1]; anns = [ann(K1)]; cls = {}
for n, (pm, ph) in enumerate(levi):
    img = rhoN(pm, ph, K1)
    hit = None
    for k, Q in enumerate(anns):
        if not F.mm(Q, img).any(): hit = k; break
    if hit is None:
        orbit.append(img); anns.append(ann(img)); hit = len(orbit) - 1
        print(f"  new kernel #{hit} from element {n}")
    cls[(pm, ph)] = hit
print(f"orbit of K_D1 under the full Levi group: {len(orbit)} kernels  ({time.time()-t0:.0f}s)")
from collections import Counter
print("  elements per kernel:", Counter(cls.values()))
# keep E_R ?
keepE = sum(1 for (pm, ph) in levi[::97] if F.rank(np.hstack([ER, rhoN(pm, ph, ER)])) != ER.shape[1])
print("  sampled Levi elements not keeping E_R:", keepE)
names = {}
for k, Kb in enumerate(orbit):
    if not F.mm(ann(K2), Kb).any(): names[k] = 'K_D2'
names[0] = 'K_D1'
for k in range(len(orbit)): names.setdefault(k, f'K_{k}')
def dimcap(idx):
    return 136 - F.rank(np.vstack([F.mm(anns[k], ER) for k in idx]))  # in E_R coords
print("  dim of each:", [orbit[k].shape[1] for k in range(len(orbit))])
for r_ in range(2, len(orbit) + 1):
    for c in itertools.combinations(range(len(orbit)), r_):
        print(f"  cap {[names[k] for k in c]}: dim {dimcap(c)}  (E_R / cap = {136 - dimcap(c)})")
pickle.dump({'cls': cls, 'names': names}, open('v_classes.pkl', 'wb'))
print(f"done {time.time()-t0:.0f}s")
