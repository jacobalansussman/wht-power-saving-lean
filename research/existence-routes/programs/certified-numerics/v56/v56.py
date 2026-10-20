# v56.py W M -- EXACT existence test for a2-gct's relaxation "W in M inside V_56" (the 56-dim irreducible of Sp(6,2)),
# using c2-finite-group's exact integral model (Ind_B^G(lambda_3) on Z^2835, signed permutations; fgcore.py, as in
# exact_cert.py) and the subspace-game reading (ex2-certnum PROOF.md Lemma 8, Corollary 11).
#  1. X = the 56 inside M (class-sum eigenvalues, exact membership), basis B of integer vectors, rows R with B_R invertible.
#  2. P_0-types: sums of A-weight spaces over GL(3,2)-orbits of characters of A (exact integer projectors 64 Pi_psi).
#     Basis T adapted to the types; every matrix rho_T(g) = T_R^-1 (rho(g) T)_R exactly (fmpq).
#  3. Checks: the P_0 generators are block diagonal in T; each type is irreducible (commutant of the generators has
#     dimension 1 modulo a prime, hence exactly 1 over Q); the types have distinct dimensions.  So the relaxation's free
#     algebra S = image of C[P_0] is ALL block-diagonal matrices (Burnside).
#  4. Exact ranks: need_tau = rank k[tau, other], cap_tau = rank h[tau, other].
#  5. Inputs-only schedule: every move's input r = sum_l beta_l e_l, chosen at random (integers) inside the linear space
#     of beta's whose deliveries into the slack types lie in those types' needed spaces; then exact rank tests that the
#     final row spaces contain every needed vector.  If they do, a W-in-M relaxed scheme with rational coefficients exists
#     (outputs by Lemma 8a/b); the outputs are then solved and the whole program checked exactly against the target.
import sys, os, time, random, json
import numpy as np
W, M = int(sys.argv[1]), int(sys.argv[2]); SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 1
TSPEC = (sys.argv[4] if len(sys.argv) > 4 else ','.join(['k'] * W)).split(',')   # per-array target: k (kernel) or h2
assert len(TSPEC) == W and all(t in ('k', 'h2') for t in TSPEC)
args = sys.argv; sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from fgcore import *
from flint import fmpq_mat, fmpz_mat, fmpq, fmpz, nmod_mat
t0 = time.time(); P31 = 2147483647
DEG = int(os.environ.get('DEG', 56)); CHI2A = int(os.environ.get('CHI2A', -24)); lam = int(os.environ.get('LAM', 3))
TAG = '%d%+d' % (DEG, CHI2A)
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
def apply(Ps, y):
    P, s = Ps; out = np.zeros_like(y); out[P] = s.astype(np.int64)[:, None] * y if y.ndim == 2 else s.astype(np.int64) * y; return out
def rank_mod(Mx, p=P31, want_pivots=False):
    Mx = np.array(Mx, dtype=np.int64) % p; r = 0; piv = []; rows, cols = Mx.shape
    for c in range(cols):
        nz = np.nonzero(Mx[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]; Mx[[r, pr]] = Mx[[pr, r]]
        iv = pow(int(Mx[r, c]), p - 2, p); Mx[r] = (Mx[r] * iv) % p
        col = Mx[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): Mx[nzr] = (Mx[nzr] - (col[nzr, None] * Mx[r][None, :]) % p) % p
        piv.append(c); r += 1
        if r == rows: break
    return (r, piv) if want_pivots else r
# ---- 1. the 56 inside M
MA = [mono(g, lam) for g in TR]; MB = [mono(g, lam) for g in C2B]
Za = np.zeros((NF, NF), np.int64); Zb = np.zeros((NF, NF), np.int64)
for P, s in MA: Za[P, np.arange(NF)] += s.astype(np.int64)
for P, s in MB: Zb[P, np.arange(NF)] += s.astype(np.int64)
alpha = 0.6180339887
w, V = np.linalg.eigh(Za + alpha * Zb)
cl = []; st = 0
for i in range(1, NF + 1):
    if i == NF or w[i] - w[i - 1] > 1e-6: cl.append((st, i)); st = i
consts = []
for (a, b) in cl:
    Q = V[:, a:b]; ra = np.trace(Q.T @ Za @ Q) / (b - a); rb = (w[a:b].mean() - ra) / alpha
    consts.append((round(ra, 6), round(rb, 6), b - a, a, b))
J = [j for j, t in enumerate(TABLE) if t[0] == DEG and t[1] == CHI2A]; assert len(J) == 1; J = J[0]
d, va, vb = TABLE[J]; z2a, z2b = 63 * va // d, 315 * vb // d
hit = [c for c in consts if abs(c[0] - z2a) < 1e-6 and abs(c[1] - z2b) < 1e-6]
assert len(hit) == 1 and hit[0][2] == d, hit
a_, b_ = hit[0][3], hit[0][4]; Q = V[:, a_:b_]
F = int(np.argmax((Q * Q).sum(axis=1))); v0 = Q @ Q[F]; c = 1451520 // d
u = np.rint(c * v0).astype(np.int64); err = np.abs(c * v0 - u).max(); assert err < 1e-6
assert (Za @ u == z2a * u).all() and (Zb @ u == z2b * u).all() and u.any()
log('irreducible %s (lambda %d): u in X checked exactly;' % (TAG, lam)); log('  (chi(2a) = %d): (z2a, z2b) = (%d, %d) unique among the constituents' % (va, z2a, z2b))
rng = random.Random(1)
def rand_elt(n=12):
    g = np.arange(NV)
    for _ in range(n): g = comp(GENS[rng.randrange(len(GENS))], g)
    return g
Bcols = [u]; rB = 1
while rB < d:
    y = apply(mono(rand_elt(), lam), u); r2 = rank_mod(np.array(Bcols + [y]))
    if r2 > rB: Bcols.append(y); rB = r2
Bm = np.array(Bcols).T
rB, R = rank_mod(Bm.T, want_pivots=True); assert rB == d
# ---- 2. types
Ael = {}
for bits in range(64):
    g = np.arange(NV)
    for i in range(6):
        if bits >> i & 1: g = comp(g, A_GENS[i])
    Ael[key(g)] = bits
def abits(g): return Ael[key(g)]
def act(l, s):
    li = inv(l); return sum(((bin(s & abits(comp(li, comp(A_GENS[j], l)))).count('1') & 1) << j) for j in range(6))
orbits = []; seen = set()
for s in range(64):
    if s in seen: continue
    orb = {s}; fr = [s]
    while fr:
        nx = []
        for x in fr:
            for l in L_GENS:
                y = act(l, x)
                if y not in orb: orb.add(y); nx.append(y)
        fr = nx
    seen |= orb; orbits.append(sorted(orb))
MAg = [mono(g, lam) for g in A_GENS]
def Pi(s, y):
    for i in range(6): y = y + (-1 if s >> i & 1 else 1) * apply(MAg[i], y)
    return y
types = []
for orb in orbits:
    cols = []; r_ = 0
    for s_ in orb:
        Y = Pi(s_, Bm)
        for cc in range(Y.shape[1]):
            yv = Y[:, cc]
            if not yv.any(): continue
            r2 = rank_mod(np.array(cols + [yv]))
            if r2 > r_: cols.append(yv); r_ = r2
    if r_: types.append((len(orb), r_, np.array(cols).T))
# ---- 3. split each orbit part into P_0-irreducibles (exactly) and check irreducibility
def fq(Mx):
    return fmpq_mat(Mx.shape[0], Mx.shape[1], [int(x) for x in Mx.ravel()])
def modp(Mq, p=P31):
    out = np.zeros((Mq.nrows(), Mq.ncols()), dtype=np.int64)
    for i in range(Mq.nrows()):
        for j in range(Mq.ncols()):
            x = Mq[i, j]; num, den = int(x.p), int(x.q); assert den % p
            out[i, j] = num % p * pow(den, p - 2, p) % p
    return out
def act_on(Tpart, g):          # matrix of rho(g) on span(Tpart) (exact), Tpart integer NF x dpart, rows Rp invertible
    dp = Tpart.shape[1]; _, Rp = rank_mod(Tpart.T, want_pivots=True)
    return fq(Tpart[Rp]).inv() * fq(apply(mono(g, lam), Tpart)[Rp])
def commutant(blocks, dt):     # exact basis (list of dt x dt fmpq_mat) of {X : X g = g X for all g}
    import math
    rows_p = []; rows_q = []
    for gq in blocks:
        gp = modp(gq)
        Cp = (np.kron(np.eye(dt, dtype=np.int64), gp.T) - np.kron(gp, np.eye(dt, dtype=np.int64))) % P31
        rows_p.append(Cp)
        # exact rows, same layout
        for a in range(dt):
            for b in range(dt):       # equation (X g - g X)[a, b] = sum_c X[a, c] g[c, b] - g[a, c] X[c, b]
                r = {}
                for c_ in range(dt):
                    if gq[c_, b] != 0: r[a * dt + c_] = r.get(a * dt + c_, 0) + gq[c_, b]
                    if gq[a, c_] != 0: r[c_ * dt + b] = r.get(c_ * dt + b, 0) - gq[a, c_]
                rows_q.append(r)
    Cp = np.concatenate(rows_p, 0)
    rk, piv = rank_mod(Cp.T, want_pivots=True)            # independent equations mod p (pivot rows of Cp)
    nul_p = dt * dt - rk
    sel = [rows_q[i] for i in piv]
    ents = []
    for r in sel:
        den = 1
        for v in r.values(): den = den * int(v.q) // math.gcd(den, int(v.q))
        row = [0] * (dt * dt)
        for k_, v in r.items(): y = v * den; row[k_] = int(y.p)
        ents.extend(row)
    Zm = fmpz_mat(len(sel), dt * dt, ents)
    X, nul = Zm.nullspace()
    assert nul == nul_p, (nul, nul_p)   # mod-p rank <= Q-rank, so nul <= nul_p; equality pins the commutant
    basis = [fmpq_mat(dt, dt, [X[i, c] for i in range(dt * dt)]) for c in range(nul)]
    return basis
pieces = []          # (dim, integer NF x dim basis)
for (osz, dimo, Tpart) in types:
    gens_blocks = [act_on(Tpart, g) for g in A_GENS + L_GENS]
    Cb = commutant(gens_blocks, dimo)
    log('orbit %d part (dim %d): commutant of P_0 has dimension %d (exact)' % (osz, dimo, len(Cb)))
    if len(Cb) == 1: pieces.append((dimo, Tpart)); continue
    assert len(Cb) == 2, 'more than two pieces: not handled'
    I_ = fmpq_mat(dimo, dimo, [1 if a == b else 0 for a in range(dimo) for b in range(dimo)])
    Zc = Cb[0] if any(Cb[0][a, b] != 0 for a in range(dimo) for b in range(dimo) if a != b) or Cb[0] != Cb[0][0, 0] * I_ else Cb[1]
    # minimal polynomial x^2 - s x - p of Zc in the 2-dim algebra {I, Zc}
    Z2 = Zc * Zc
    # solve Z2 = s Zc + p I using two entries
    import itertools as it_
    found = None
    for (a1, b1), (a2, b2) in it_.combinations([(a, b) for a in range(dimo) for b in range(dimo)], 2):
        m = fmpq_mat([[Zc[a1, b1], I_[a1, b1]], [Zc[a2, b2], I_[a2, b2]]])
        if m.det() != 0:
            sol = m.solve(fmpq_mat([[Z2[a1, b1]], [Z2[a2, b2]]])); found = (sol[0, 0], sol[1, 0]); break
    s_, p_ = found
    assert Z2 == s_ * Zc + p_ * I_, 'not quadratic'
    disc = s_ * s_ + 4 * p_
    num, den = int(disc.p), int(disc.q)
    import math
    rn, rd = math.isqrt(num), math.isqrt(den); assert rn * rn == num and rd * rd == den, 'eigenvalues not rational'
    l1 = (s_ + fmpq(rn, rd)) / 2; l2 = (s_ - fmpq(rn, rd)) / 2
    E1 = (Zc - l2 * I_) / (l1 - l2); E2 = I_ - E1
    assert E1 * E1 == E1 and E2 * E2 == E2
    for E in (E1, E2):
        r_ = E.rank(); Ecols = []
        # independent columns of E (exact), scaled to integers, mapped back to NF-dim vectors
        sel = []; 
        for cidx in range(dimo):
            cand = sel + [cidx]
            sub_ = fmpq_mat(dimo, len(cand), [E[a, c_] for a in range(dimo) for c_ in cand])
            if sub_.rank() == len(cand): sel = cand
            if len(sel) == r_: break
        cols_int = []
        for c_ in sel:
            den = 1
            for a in range(dimo): den = den * int(E[a, c_].q) // math.gcd(den, int(E[a, c_].q))
            v = np.array([int((E[a, c_] * den).p) for a in range(dimo)], dtype=object)
            cols_int.append(Tpart.astype(object).dot(v))
        pieces.append((r_, np.array(cols_int, dtype=object).T))
        log('  split: piece of dim %d' % r_)
pieces.sort(key=lambda t: t[0])
dims = [p[0] for p in pieces]; off = np.cumsum([0] + dims)
log('P_0-irreducible pieces of %s:' % TAG, dims)
assert sum(dims) == d and len(set(dims)) == len(dims), 'pieces must have distinct dimensions (then pairwise non-isomorphic)'
T = np.concatenate([p[1] for p in pieces], 1)
Tm = np.array([[int(x) % P31 for x in row] for row in T], dtype=np.int64)
_, R = rank_mod(Tm.T, want_pivots=True); assert len(R) == d
def fqo(Mx): return fmpq_mat(Mx.shape[0], Mx.shape[1], [int(x) for x in np.asarray(Mx).ravel()])
TRi = fqo(T[R]).inv()
def rhoT(g):
    P, s = mono(g, lam); out = np.empty_like(T); out[P] = (s.astype(np.int64).astype(object))[:, None] * T
    return TRi * fqo(out[R])
gensP0 = [rhoT(g) for g in A_GENS + L_GENS]
def blockdiag(Mq):
    for i in range(d):
        for j in range(d):
            bi = np.searchsorted(off, i, side='right') - 1; bj = np.searchsorted(off, j, side='right') - 1
            if bi != bj and Mq[i, j] != 0: return False
    return True
assert all(blockdiag(g) for g in gensP0), 'a P_0 generator is not block diagonal'
for ti, dt in enumerate(dims):
    blocks = [fmpq_mat(dt, dt, [g[a, b] for a in range(off[ti], off[ti + 1]) for b in range(off[ti], off[ti + 1])]) for g in gensP0]
    nc = len(commutant(blocks, dt))
    log('piece %d (dim %d): commutant of the P_0 generators has dimension %d (exact)' % (ti, dt, nc))
    assert nc == 1
h = rhoT(h1); k = rhoT(k3); hh2 = rhoT(h2); TG = {'k': k, 'h2': hh2}
def blk(Mq, a, b):   # rows of type a, columns of the other types (b = list of types)
    rows = list(range(off[a], off[a + 1])); cols = [j for t in b for j in range(off[t], off[t + 1])]
    return fmpq_mat(len(rows), len(cols), [Mq[i, j] for i in rows for j in cols])
# dump the exact data (adapted basis) for later analysis, and the block ranks of h, h2, k
json.dump({'dims': dims, 'h': [[str(h[i, j]) for j in range(d)] for i in range(d)], 'k': [[str(k[i, j]) for j in range(d)] for i in range(d)],
           'h2': [[str(hh2[i, j]) for j in range(d)] for i in range(d)]}, open('exact_data_%s.json' % TAG, 'w'))
for nm_, Mq_ in [('h', h), ('k', k), ('h2', hh2)]:
    rk_ = [[fmpq_mat(dims[a_], dims[b_], [Mq_[i, j] for i in range(off[a_], off[a_ + 1]) for j in range(off[b_], off[b_ + 1])]).rank()
            for b_ in range(len(dims))] for a_ in range(len(dims))]
    log('block ranks of %s (rows piece a, columns piece b):' % nm_, rk_)
nt = len(dims); needk = {}; cap = []
for a in range(nt):
    oth = [t for t in range(nt) if t != a]
    cap.append(blk(h, a, oth).rank())
    for nm, Tg in TG.items(): needk[(nm, a)] = blk(Tg, a, oth).rank()
need = [sum(needk[(TSPEC[l], a)] for l in range(W)) for a in range(nt)]          # total over the arrays
log('exact: need per array (k)', [needk[('k', a)] for a in range(nt)], ' (h2)', [needk[('h2', a)] for a in range(nt)],
    ' cap per move', cap, '; targets', TSPEC, '; total need', need, ' moves needed >=', [round(n / c, 4) for n, c in zip(need, cap)])
for a in range(nt):
    if M * cap[a] < need[a]:
        log('CAPACITY BOUND VIOLATED for type %d: %d moves x %d < total need %d -> NO relaxed scheme of this shape' % (a, M, cap[a], need[a])); sys.exit()
# ---- 5. inputs-only schedule.  A move's input r = sum_l beta_l e_l with beta_l = diag(c_0, ..., c_{nt-1}) (block
# diagonal, so in S).  Into type a it delivers the rows of sum_{b != a} h[a,b] c_b^(l) in the type-b coordinates of array l.
def full(Mq): return np.array([[Mq[i, j] for j in range(Mq.ncols())] for i in range(Mq.nrows())], dtype=object)
hF = full(h); TF = {nm: full(Tg) for nm, Tg in TG.items()}
def sub(Ff, a, b): return Ff[off[a]:off[a + 1], off[b]:off[b + 1]]
aim = {(nm, a) for nm in TG for a in range(nt) if needk[(nm, a)] < sum(dims) - dims[a]}   # proper needed subspaces
# needed space of type a (per array): rowspace of k[a, others] in C^{others}; its annihilator Ann_a (columns)
def nullspace_q(Mq):     # right null space basis of an fmpq_mat, as an object matrix (columns)
    import math
    den = 1
    for i in range(Mq.nrows()):
        for j in range(Mq.ncols()):
            q_ = int(Mq[i, j].q); den = den * q_ // math.gcd(den, q_)
    ents = []
    for i in range(Mq.nrows()):
        for j in range(Mq.ncols()):
            y = Mq[i, j] * den; assert int(y.q) == 1; ents.append(int(y.p))
    Z = fmpz_mat(Mq.nrows(), Mq.ncols(), ents)
    X, nul = Z.nullspace()
    return np.array([[int(X[i, j]) for j in range(nul)] for i in range(X.nrows())], dtype=object)
ann = {}
for (nm, a) in aim:
    oth = [t for t in range(nt) if t != a]; ann[(nm, a)] = nullspace_q(blk(TG[nm], a, oth))
log('needed subspaces to aim at (target, type):', sorted(aim), ' annihilator dims', {key_: v.shape[1] for key_, v in ann.items()})
# unknowns x = (vec c_0, vec c_1, ...) (row-major), constraints: for each slack type a, [h[a,b] c_b]_{b != a} Ann_a = 0
nunk = sum(dt * dt for dt in dims); uoff = np.cumsum([0] + [dt * dt for dt in dims])
AIMALL = os.environ.get('AIMALL') == '1'
zr = random.Random(SEED + 1000)
ZS = {}
for a in range(nt):      # t_a aimed left vectors for piece a: enough for the total need, at most the capacity
    t_a = min(cap[a], -(-need[a] // M) + int(os.environ.get('AIMEXTRA', '1')))
    ZS[a] = [np.array([zr.randint(-3, 3) for _ in range(dims[a])], dtype=object) for _ in range(t_a)]
log('aimed left vectors per move:', {a: len(ZS[a]) for a in range(nt)}, '(AIMALL)' if AIMALL else '')
def constraint_basis(nm):
  eqs = []
  for a in [a_ for (n_, a_) in aim if n_ == nm]:
    oth = [t for t in range(nt) if t != a]
    Ann = ann[(nm, a)]; coff = np.cumsum([0] + [dims[t] for t in oth])
    # AIM (junk-tolerant): only t_a chosen left vectors z per move must deliver into the needed space; the other
    # directions of the move deliver junk, which is harmless (containment).  AIMALL=1 restores 'every row aimed'.
    zlist = [np.eye(dims[a], dtype=object)[i] for i in range(dims[a])] if AIMALL else ZS[a]
    for z in zlist:
        for q in range(Ann.shape[1]):
            row = [0] * nunk
            for bi, b in enumerate(oth):
                zH = z.dot(sub(hF, a, b)); Ab = Ann[coff[bi]:coff[bi + 1], q]
                # sum_{u, v} (z H)[u] c_b[u, v] Ann[v, q]
                for uu in range(dims[b]):
                    if zH[uu] == 0: continue
                    for vv in range(dims[b]):
                        if Ab[vv] != 0: row[uoff[b] + uu * dims[b] + vv] += zH[uu] * Ab[vv]
            eqs.append(row)
  log('constraint system for %s-arrays: %d x %d' % (nm, len(eqs), nunk))
  if eqs:
    Eq = fmpq_mat(len(eqs), nunk, [fmpq(x) if not isinstance(x, fmpq) else x for r in eqs for x in r])
    NB = nullspace_q(Eq)
  else:
    NB = np.eye(nunk, dtype=object)
  log('free parameters per (move, %s-array): %d' % (nm, NB.shape[1]))
  return NB
NBS = {nm: constraint_basis(nm) for nm in set(TSPEC)}
rr = random.Random(SEED)
def draw(nm):
    NB = NBS[nm]
    coef = [rr.randint(-3, 3) for _ in range(NB.shape[1])]
    x = NB.dot(np.array(coef, dtype=object))
    return [np.array(x[uoff[t]:uoff[t + 1]], dtype=object).reshape(dims[t], dims[t]) for t in range(nt)]
moves = [[draw(TSPEC[l]) for l in range(W)] for j in range(M)]
# deliveries and final containment, per type, exactly
ok = True
for a in range(nt):
    oth = [t for t in range(nt) if t != a]
    width = sum(dims[t] for t in oth)
    Drows = []
    for j in range(M):
        Dj = np.zeros((dims[a], width * W), dtype=object)
        for l in range(W):
            co = 0
            for b in oth:
                Dj[:, l * width + co: l * width + co + dims[b]] = sub(hF, a, b).dot(moves[j][l][b]); co += dims[b]
        Drows.append(Dj)
    D = np.concatenate(Drows, 0)
    Nd = np.zeros((dims[a] * W, width * W), dtype=object)
    for l in range(W):
        co = 0
        for b in oth:
            Nd[l * dims[a]:(l + 1) * dims[a], l * width + co: l * width + co + dims[b]] = sub(TF[TSPEC[l]], a, b); co += dims[b]
    def toq(Ar): return fmpq_mat(Ar.shape[0], Ar.shape[1], [fmpq(x) if not isinstance(x, fmpq) else x for x in Ar.ravel()])
    rD = toq(D).rank(); rDN = toq(np.concatenate([D, Nd], 0)).rank()
    log('type %d: delivered rank %d (capacity %d x %d = %d, total need %d); with the needed rows %d -> %s' %
        (a, rD, M, cap[a], M * cap[a], need[a], rDN, 'CONTAINS the needed space' if rD == rDN else 'missing %d' % (rDN - rD)))
    ok &= rD == rDN
print('VERDICT: %s' % (('EXISTS: a %d-in-%d relaxed scheme in the irreducible %s with rational coefficients (inputs-only moves; the outputs exist by '
      'the module lemma since every needed row lies in the row spaces)' % (W, M, TAG)) if ok else 'not closed by this draw'), flush=True)
json.dump({'W': W, 'M': M, 'seed': SEED, 'targets': TSPEC, 'dims': dims, 'need': need, 'cap': cap, 'ok': bool(ok)}, open('out_%s_%din%d_%s_s%d.json' % (TAG, W, M, '-'.join(TSPEC), SEED), 'w'))
