"""vd4addr.py -- v-e3-bpexact (skeptic): the same own rebuild as vd4V.py, in the ADDRESS-TYPE representations gl15 (linear
action on the 15 nonzero addresses = 1 + gl14) and aff16 (affine permutation of the 16 addresses), for the generic lift.
slot-c coefficient of a product M = g_d p_j k_b f0 is beta = M rep_c^-1 (rep_c: a fixed element of slot c); sigma(beta) from
beta's address permutation (own code).
(Header of vd4V.py follows.)

In V the image of a cell's formal sum is just the ACTUAL matrix sum of its terms, so the V test needs no slot-module
machinery: the top's F1 part in character a is  X_a . D3 . A . P_a . kappa . f0.  Upper-bound ('fullrow') version: the
top's coefficient is replaced by an arbitrary row block W = (W_d)_d over the 30 cosets g_d = D3 b_d of B D3 B (my own
coset reps, from my own letter sample); the slot-c sum is  sum_{d,j,b: cell(g_d p_j k_b f0) = c} chi_a(j)/7 y_b W_d (g_d p_j k_b f0).
Slots of products are computed DIRECTLY from the product matrices (no cocycles, no permutation composition).
Also: legitimacy of every letter (Z4-monomial, Clifford, fixes L_Z), m = D_z * letter, absorption by my own formal-sum code,
the flat lift (value g_d p_j f0, the k_b's Pauli correction), and a second top direction D3.
argv: seed psi nBm D3z"""
import sys, os, time, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lift import *          # Lift, mm, minv, key, I16, G, sigma, LZ, img_mask, table, ...
t0 = time.time()
seed, psi, nBm, z3 = (int(x) for x in (sys.argv[1:] + ['21', '1', '400', '1'][len(sys.argv) - 1:])[:4])
L = Lift(seed=seed, psi=psi, verbose=True)
P = L.P
# --- the same letter list as run1.py / fullrow.py (needed to match the stored designed ys)
rng = random.Random(seed + 1)
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < nBm:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
ys = np.load(os.path.join(HERE, f'ys_seed{seed}_psi{psi}.npy'))
assert len(ys) == len(Bm), (len(ys), len(Bm))
mi = minv(L.m); ks = [mm(L.m, b, mi) for b in Bm]
if os.environ.get('CTL'):      # POSITIVE CONTROL: the plain chain (kappa = 1, no twist): the top's D3 move on f0 reaches K
    Bm = [I16]; ys = np.array([1.0 + 0j]); ks = [I16]; L.p = [I16] * 7; L.qp = [I16] * 7
    print('[CTL] plain-chain positive control: kappa = 1, p_j = 1 (only a = 0 is nonzero)', flush=True)

# ---------------------------------------------------------------- (1) legitimacy, my own checks
def z4mono(M):
    nz = np.abs(M) > 1e-9
    if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()): return False
    v = M[nz]
    return bool(np.all(np.min(np.abs(v[:, None] - np.array([1, 1j, -1, -1j])[None, :]), axis=1) < 1e-9))
def fixesLZ(M):
    return img_mask(table(sigma(M)), LZ) == LZ          # sigma asserts Clifford
bad = []
named = [('c', L.c), ("c'", L.cp), ('p', L.p), ('beta', L.beta), ('q', [L.q])]
for nm, lst in named:
    for i, X in enumerate(lst):
        if not (z4mono(X) and fixesLZ(X)): bad.append((nm, i))
# the coefficient letters F1 actually applies to O2_h:  q^h b c_h^-1  (F1 = m . sum psi^-1 y_b (q^h b c_h^-1) O2_h)
nco = 0
for h in range(7):
    for b in Bm:
        X = mm(L.qp[h], b, minv(L.c[h])); nco += 1
        if not (z4mono(X) and fixesLZ(X)): bad.append(('coef', h))
# b in B_m fixes Lam_m and the k_b = m b m^-1 are letters
for b, k in zip(Bm, ks):
    if not (z4mono(k) and fixesLZ(k)): bad.append(('k', 0))
# input symmetry: c_h c'_h must be a letter (tau_h maps inputs to letter * input)
for h in range(7):
    if not (z4mono(mm(L.c[h], L.cp[h])) and fixesLZ(mm(L.c[h], L.cp[h]))): bad.append(('cc', h))
# m = D_z * letter for some z; m1, m2 likewise
def is_move(M):
    for z in range(1, 16):
        X = mm(minv(G.D[z]), M)
        if z4mono(X) and fixesLZ(X): return z
    return None
print(f'[legit] letters checked: c, c\', p, beta, q, {nco} F1 coefficient letters, {len(ks)} k_b, c_h c\'_h; bad = {bad[:5]} ({len(bad)})')
print(f'[legit] m = D_{is_move(L.m)} * letter; m1 = D_{is_move(L.m1)} * letter; m2 = D_{is_move(L.m2)} * letter;'
      f' c^7 = I {np.allclose(np.linalg.matrix_power(L.c[1], 7), I16)}; q^7 = I {np.allclose(np.linalg.matrix_power(L.q, 7), I16)}')
print(f'[legit] p homomorphism: {all(np.allclose(mm(L.p[a], L.p[b]), L.p[(a + b) % 7]) for a in range(7) for b in range(7))}')

# ---------------------------------------------------------------- (2) absorption, my own formal sums
def fsum_add(u, l, w, c):
    a = np.round(w * 4096); k = (l % 7, np.concatenate([a.real, a.imag]).astype(np.int64).tobytes())
    if k in u: u[k][0] += c
    else: u[k] = [c, w]
F1 = {}
for h in range(7):
    ps = np.exp(-2j * np.pi * psi * h / 7)
    for y, b in zip(ys, Bm):
        X = mm(L.qp[h], b, minv(L.c[h]))                 # coefficient letter
        O2h = mm(L.c[h], L.g2, L.cp[h])                   # tau_h(O2) = [c_h g2 c'_h] x_h
        fsum_add(F1, h, mm(L.m, X, O2h), ps * y)
worst = 0.0
for g in range(1, 7):
    lhs = {}; rhs = {}
    for (l, kk), (c, w) in F1.items():
        fsum_add(lhs, l + g, mm(L.c[g], w, L.cp[g]), c)
        fsum_add(rhs, l, mm(L.c[g], minv(L.p[g]), w), np.exp(2j * np.pi * psi * g / 7) * c)
    keys = set(lhs) | set(rhs)
    worst = max(worst, max(abs((lhs[k][0] if k in lhs else 0) - (rhs[k][0] if k in rhs else 0)) for k in keys))
sc = max(abs(v[0]) for v in F1.values())
print(f'[absorb] F1 terms {len(F1)}; max |tau_g(F1) - psi(g) c_g p_g^-1 F1| = {worst:.1e} (coef scale {sc:.2f})')
cells0 = set(L.cell(v[1]) for k, v in F1.items() if k[0] == 0)
print(f'[absorb] cells of array-0 terms {sorted(cells0)} (F = {L.F}); kappa f0 check:',
      max(abs(v[0]) for v in F1.values()) > 0, f'{time.time() - t0:.0f}s', flush=True)

# ---------------------------------------------------------------- (3) the V test, own construction
D3 = G.D[z3]
N3 = L.cell(D3); nbZ = list(P.nb[L.iZ]); assert N3 in nbZ
# own coset reps: g_d = D3 b over my own letters (fresh sample, own rng)
myrng = random.Random(1000 + seed + z3)
reps = {}
reps[N3] = D3
while len(reps) < 30:
    b = G.rand_engine(myrng); g = mm(D3, b); c = L.cell(g)
    if c not in reps: reps[c] = g
cos = sorted(reps); gd = [reps[c] for c in cos]
nbF = sorted(P.nb[L.F])
f0 = L.f0
slot_of = {}            # slot index by the CELL of the product (a neighbour of F)
for i, c in enumerate(nbF): slot_of[c] = i
iK = L.iK; assert iK in slot_of
d0 = slot_of[iK]
S4 = [slot_of[c] for c in nbF if P.dZ[c] == 4]
Sx = [s for s in S4 if s != d0]
print(f'[V] D3 = D_{z3}; |S4| {len(S4)}; d0 slot {d0}', flush=True)
Kf = np.array([mm(k, f0) for k in ks])                       # 400 x 16 x 16
Kc = [L.cell(mm(k, f0)) for k in ks]
assert all(c == L.F for c in Kc)
# symplectic images (sigma of each factor; sigma asserts Clifford); cell(w) = sigma(w)^-1 L_Z, sigma multiplicative
Sf0 = sigma(f0); Sk = [mul(sigma(k), Sf0) for k in ks]
Sgp = [[sigma(mm(g, L.p[j])) for j in range(7)] for g in gd]
SkI = [inv(x) for x in Sk]; SgpI = [[inv(x) for x in row] for row in Sgp]
_LZb = [1 << (4 + j) for j in range(4)]
def cellS(A, B, AI=None, BI=None):        # cell of the product with symplectic image A B:  (A B)^-1 L_Z = B^-1 A^-1 L_Z
    AI = inv(A) if AI is None else AI; BI = inv(B) if BI is None else BI
    return P.idx[span_mask([app(BI, app(AI, v)) for v in _LZb])]
# spot check: the symplectic composition gives the same cell as the product matrix
for t in range(20):
    d, j, ib = myrng.randrange(30), myrng.randrange(7), myrng.randrange(len(ks))
    assert cellS(Sgp[d][j], Sk[ib]) == L.cell(mm(gd[d], L.p[j], ks[ib], f0))
chars = [[np.exp(2j * np.pi * a * j / 7) / 7 for j in range(7)] for a in range(7)]
SL = [[[cellS(None, None, SgpI[d][j], SkI[ib]) for ib in range(len(ks))] for j in range(7)] for d in range(30)]
print(f'[V] slots of all 30 x 7 x {len(ks)} products  {time.time() - t0:.0f}s', flush=True)

# slot representatives: the first product found in each slot
rep = {}
for d in range(30):
    for j in range(7):
        for ib in range(len(ks)):
            c = slot_of[SL[d][j][ib]]
            if c not in rep: rep[c] = mm(gd[d], L.p[j], ks[ib], f0)
assert len(rep) == 30
repi = {c: minv(r) for c, r in rep.items()}
def addr_perm(Bt):          # batch of monomial matrices -> image index of each address y (column y has its 1 in row pi(y))
    return np.argmax(np.abs(Bt), axis=1)
def sig_of(name, pi):
    if name == 'gl15':
        M = np.zeros((15, 15)); c0 = pi[0]
        for v in range(1, 16): M[(pi[v] ^ c0) - 1, v - 1] = 1
        return M
    if name == 'aff16':
        M = np.zeros((16, 16))
        for y in range(16): M[pi[y], y] = 1
        return M
    if name == 'triv': return np.ones((1, 1))
def buildS(name, ds, yv=None):
    yv = ys if yv is None else yv
    Ma = {a: np.zeros((30 * ds, 30 * ds), dtype=complex) for a in range(7)}
    cache = {}
    for d in range(30):
        for j in range(7):
            Gdj = gd[d] @ L.p[j]
            vals = np.einsum('ik,bkl->bil', Gdj, Kf)
            for ib in range(len(ks)):
                c = slot_of[SL[d][j][ib]]
                beta = vals[ib] @ repi[c]
                pi = tuple(int(x) for x in np.argmax(np.abs(beta), axis=0))
                assert z4mono(np.round(beta * 4096) / 4096)
                if pi not in cache: cache[pi] = sig_of(name, pi)
                v = cache[pi]
                for a in range(7):
                    Ma[a][d * ds:(d + 1) * ds, c * ds:(c + 1) * ds] += chars[a][j] * yv[ib] * v
    return Ma
def ranks(Ma, cons, label, ds, mlist):
    colx = np.concatenate([np.arange(s * ds, (s + 1) * ds) for s in cons]); col0 = np.arange(d0 * ds, (d0 + 1) * ds)
    out = []; scale = max(np.abs(Ma[a]).max() for a in range(7))
    for a in range(7):
        C = Ma[a][:, colx]; U, s, Vh = np.linalg.svd(C, full_matrices=True); r = int((s > 1e-9 * scale).sum())
        Nul = U[:, r:].conj().T; junk = float(np.abs(Nul @ C).max()) if Nul.shape[0] else 0.0
        B0 = Nul @ Ma[a][:, col0]; sv = np.linalg.svd(B0, compute_uv=False) if B0.size else np.zeros(ds)
        ra = int((sv > 1e-8 * scale).sum()); gap = (sv[ra - 1] if ra else 0, sv[ra] if ra < len(sv) else 0)
        out.append((a, ra, mlist[a], junk, gap))
    print(f'[{label}] ' + '; '.join(f'a{a}: r={ra} m={ma} junk={jk:.0e} sv[{g[0]:.1e}|{g[1]:.1e}]' for a, ra, ma, jk, g in out), flush=True)
    print(f'[{label}] total deficit {sum(ds - x[1] for x in out)} vs ds {ds}', flush=True)
    return out
for name, ds in (('triv', 1), ('gl15', 15), ('aff16', 16)):
    # m_a(p) in sigma: multiplicity of eigenvalue for p's own sigma-image
    sp = []
    for j in range(7):
        pi = tuple(int(x) for x in np.argmax(np.abs(L.p[j]), axis=0)); sp.append(sig_of(name, pi))
    ml = [int(round(np.trace(sum(chars[a][j] * sp[j] for j in range(7))).real)) for a in range(7)]
    Ma = buildS(name, ds)
    ranks(Ma, Sx, f'{name} generic d4', ds, ml)
print(f'done {time.time() - t0:.0f}s')
