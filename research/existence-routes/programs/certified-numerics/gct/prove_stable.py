# prove_stable.py (ex2-certnum): Proposition 6(a) in action -- a numerical-only existence proof at a STABLE point.
# Object: a2-gct's single-irreducible relaxation, W = 1 array, M paid moves, irreducible IRR (types tl from relax2.load).
# Polynomial form (no normalisation, no projection; exact in the data h, k):
#   r_1 = b_1,  u_1 = h r_1,  r_i = b_i + sum_{j<i} g_ij u_j,  u_i = h r_i,
#   R  =  k - Z - sum_i Y_i u_i  = 0          (n x n equations)
# with Z, Y_i, b_i, g_ij in S (block-diagonal: one d x d block per type, repeated on its copies).
# Steps: (1) numerical solution (Gauss-Newton, min-norm steps); (2) Jacobian rank = n^2 (stability: full row rank);
# (3) pivoted columns -> square n^2 x n^2 system, other unknowns fixed at their float values (exact binary numbers);
# (4) Krawczyk test in Arb with the DATA h, k as balls of radius DELTA around the stored floats.  If it passes:
#   for every h', k' within DELTA (entrywise) of the stored matrices, an exact relaxed scheme exists.  [no exact closure]
import sys, time, numpy as np
import relax2
from flint import arb, arb_mat, ctx
IRR, M = int(sys.argv[1]), int(sys.argv[2]); DELTA = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-12
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 1
n, h, k, tl = relax2.load([IRR], np.random.default_rng(12345))
# S coordinates: (type, a, b) -> list of matrix positions (one per copy)
coords = []
for d, slots in tl:
    for a in range(d):
        for b in range(d):
            coords.append([(s[a], s[b]) for s in slots])
sd = len(coords)
def smat(p):
    S = np.zeros((n, n))
    for c, pos in zip(p, coords):
        for (i, j) in pos: S[i, j] = c
    return S
pairs = [(i, j) for i in range(M) for j in range(i)]
names = ['Z'] + ['Y%d' % i for i in range(M)] + ['b%d' % i for i in range(M)] + ['g%d%d' % p for p in pairs]
nblk = len(names); NU = nblk * sd
def unpack(u):
    return {nm: smat(u[q * sd:(q + 1) * sd]) for q, nm in enumerate(names)}
def forward(B, hh):
    U = []; Rr = []
    for i in range(M):
        r = B['b%d' % i].copy()
        for j in range(i): r = r + B['g%d%d' % (i, j)] @ U[j]
        Rr.append(r); U.append(hh @ r)
    return U, Rr
def resid(u):
    B = unpack(u); U, _ = forward(B, h)
    R = k - B['Z'] - sum(B['Y%d' % i] @ U[i] for i in range(M))
    return R.ravel()
def jac(u):
    # numerical Jacobian by forward-mode on the S coordinates (exact formulas)
    B = unpack(u); U, Rr = forward(B, h)
    J = np.zeros((n * n, NU))
    # dU_i / d(param) propagated: build for each block parameter the derivative of every U
    for q, nm in enumerate(names):
        for c, pos in enumerate(coords):
            E = np.zeros((n, n))
            for (i, j) in pos: E[i, j] = 1.0
            if nm == 'Z': dR = -E
            elif nm[0] == 'Y': dR = -E @ U[int(nm[1:])]
            else:
                dU = []
                for i in range(M):
                    dr = np.zeros((n, n))
                    if nm == 'b%d' % i: dr = dr + E
                    for j in range(i):
                        if nm == 'g%d%d' % (i, j): dr = dr + E @ U[j]
                        dr = dr + B['g%d%d' % (i, j)] @ dU[j]
                    dU.append(h @ dr)
                dR = -sum(B['Y%d' % i] @ dU[i] for i in range(M))
            J[:, q * sd + c] = dR.ravel()
    return J
rng = np.random.default_rng(SEED)
t0 = time.time()
best = None
for st in range(6):
    u = rng.standard_normal(NU) * 0.5
    for it in range(60):
        r = resid(u); nr = np.linalg.norm(r)
        if nr < 1e-13: break
        J = jac(u); du = np.linalg.lstsq(J, -r, rcond=None)[0]
        lam = 1.0
        while lam > 1e-4 and np.linalg.norm(resid(u + lam * du)) > nr: lam /= 2
        u = u + lam * du
    if best is None or nr < best[1]: best = (u.copy(), nr)
    if nr < 1e-13: break
u, nr = best
J = jac(u); sv = np.linalg.svd(J, compute_uv=False); rk = int((sv > 1e-9 * sv[0]).sum())
print('irr %d n=%d, 1 in %d: %d unknowns, %d equations; Gauss-Newton residual %.2e; Jacobian rank %d (sv %.2e .. %.2e) [%.0fs]'
      % (IRR, n, M, NU, n * n, nr, rk, sv[0], sv[min(len(sv), n * n) - 1], time.time() - t0), flush=True)
if nr > 1e-12 or rk < n * n:
    print('NOT STABLE or not solved: no proof attempted'); sys.exit()
# square slice: pivoted QR on J's columns
import scipy.linalg as sla
_, _, piv = sla.qr(J, pivoting=True, mode='economic')
cols = sorted(piv[:n * n]); Js = J[:, cols]
print('square slice: cond %.2e' % np.linalg.cond(Js), flush=True)
# Newton on the square system (other unknowns fixed at their float values)
for it in range(5):
    r = resid(u); du = np.linalg.solve(J[:, cols], -r); u[cols] += du; J = jac(u)
print('after square Newton: residual %.2e' % np.linalg.norm(resid(u)), flush=True)
# ---------------- Krawczyk in Arb with data balls
ctx.prec = 128
def ball(x, r=0.0): return arb(float(x), r) if r else arb(float(x))
hB = arb_mat([[ball(h[i, j], DELTA) for j in range(n)] for i in range(n)])
kB = arb_mat([[ball(k[i, j], DELTA) for j in range(n)] for i in range(n)])
def smatB(vals):     # vals: list of arb, length sd
    S = arb_mat(n, n)
    for c, pos in zip(vals, coords):
        for (i, j) in pos: S[i, j] = c
    return S
def unpackB(vals):
    return {nm: smatB(vals[q * sd:(q + 1) * sd]) for q, nm in enumerate(names)}
def forwardB(B):
    U = []
    for i in range(M):
        r = B['b%d' % i]
        for j in range(i): r = r + B['g%d%d' % (i, j)] * U[j]
        U.append(hB * r)
    return U
def residB(vals):
    B = unpackB(vals); U = forwardB(B)
    R = kB - B['Z']
    for i in range(M): R = R - B['Y%d' % i] * U[i]
    return [R[i, j] for i in range(n) for j in range(n)], B, U
def jacB(vals, B, U):
    # columns for the chosen unknowns only, with all quantities as balls over the box
    out = arb_mat(n * n, len(cols))
    for ci, gidx in enumerate(cols):
        q, c = divmod(gidx, sd); nm = names[q]; pos = coords[c]
        E = arb_mat(n, n)
        for (i, j) in pos: E[i, j] = 1
        if nm == 'Z': dR = -E
        elif nm[0] == 'Y': dR = -(E * U[int(nm[1:])])
        else:
            dU = []
            for i in range(M):
                dr = arb_mat(n, n)
                if nm == 'b%d' % i: dr = dr + E
                for j in range(i):
                    if nm == 'g%d%d' % (i, j): dr = dr + E * U[j]
                    dr = dr + B['g%d%d' % (i, j)] * dU[j]
                dU.append(hB * dr)
            dR = arb_mat(n, n)
            for i in range(M): dR = dR - B['Y%d' % i] * dU[i]
        for a in range(n):
            for b in range(n): out[a * n + b, ci] = dR[a, b]
    return out
y = [arb(float(v)) for v in u]
Ry, _, _ = residB(y)
Cn = np.linalg.inv(Js)
C = arb_mat([[arb(float(Cn[i, j])) for j in range(n * n)] for i in range(n * n)])
t1 = time.time()
for rho in [1e-11, 1e-9, 1e-7]:
    X = list(y)
    for gidx in cols: X[gidx] = arb(float(u[gidx]), rho)
    RX, BX, UX = residB(X)
    JX = jacB(X, BX, UX)
    CJ = C * JX
    I = arb_mat([[1 if i == j else 0 for j in range(n * n)] for i in range(n * n)])
    D = arb_mat([[arb(0, rho)] for _ in range(n * n)])
    Ky = C * arb_mat([[v] for v in Ry])
    K = [-Ky[i, 0] for i in range(n * n)]
    T2 = (I - CJ) * D
    ok = True; worst = 0.0
    for i in range(n * n):
        kv = K[i] + T2[i, 0]          # displacement of K(X) from y, must lie strictly inside [-rho, rho]
        if not arb(0, rho).contains_interior(kv): ok = False
        worst = max(worst, float(abs(kv).upper()))
    print('Krawczyk at half-width %.0e: %s (max |K - y| %.2e vs %.0e) [%.0fs]' % (rho, 'PASS' if ok else 'fail', worst, rho, time.time() - t1), flush=True)
    if ok:
        print('THEOREM [certified numerics, no exact closure]: for every h\', k\' within %.0e of the stored irr-%d data, the'
              ' relaxed 1-in-%d scheme exists exactly (unique in the box for the %d chosen unknowns, the others fixed).'
              % (DELTA, IRR, M, n * n), flush=True)
        break
