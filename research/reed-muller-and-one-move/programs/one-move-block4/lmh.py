#!/usr/bin/env python3
"""lmh.py (t2-block4): complex Levenberg-Marquardt with an analytic Jacobian for the NO-PHASE engine game with a finite
address group H, any H (not only cyclic).  Exact formulation by e3-rmbuild Theorem A + Lemma C:
  a W-in-M unit exists iff for every binary quadratic form Q (one per H-orbit suffices)
      Atil (I - Pi_Q Gamma)^-1 Pi_Q Btil = I_W (x) O_Q                                              (*)
  in V_Q = C^H (channels k in H).  Coefficients are elements of C[H] acting LEFT-regularly (P_g |k> = |gk>);
  Pi_Q = diag over (move t, channel k) of [Q(k^-1 z_t) = 1]; O_Q = diag over k of [tr(Q o k^-1) = 1].
  Derivation (mine, same as Lemma C): y_t = Pi_t (B_t + sum_{s<t} Gamma_ts y_s) is the new information of move t
  ([zz^T] = 1 - 2 Pi in V_Q); the trivial form forces the free part of the output to be the input.
Unknowns: A (W x M x n), B (M x W x n), G (strictly lower M x M x n), complex.
Residual is holomorphic, so complex Gauss-Newton:  (J^H J + lam D) dx = -J^H r.
Usage: lmh.py --H <spec> --W 3 --M 11 [--dirs z1,z2,...] [--plain] [--restarts R] [--seed s] [--check]
H spec: trivial | singer:k | gens:r0.r1.r2.r3/r0.r1.r2.r3 (row bitmasks of generators, '/' between generators)
"""
import argparse, itertools, time, json, sys
import numpy as np
np.seterr(all="ignore")

m = 4
def parity(x): return bin(x).count('1') & 1
def apply(h, x): return sum(parity(h[i] & x) << i for i in range(m))
def matmul(A, B):
    out = []
    for i in range(m):
        r = 0
        for k in range(m):
            if (A[i] >> k) & 1: r ^= B[k]
        out.append(r)
    return tuple(out)
IDM = tuple(1 << i for i in range(m))
IDX = [(i, j) for i in range(m) for j in range(i, m)]
NF = 1 << len(IDX)
DIRS = list(range(1, 1 << m))
def Qval(Q, x):
    v = 0
    for b, (i, j) in enumerate(IDX):
        if (Q >> b) & 1 and (x >> i) & 1 and (x >> j) & 1: v ^= 1
    return v
QV = np.array([[Qval(Q, x) if x else 0 for x in range(1 << m)] for Q in range(NF)], dtype=np.int8)

def gen_group(gens):
    grp = {IDM}; frontier = [IDM]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = matmul(s, g)
                if h not in grp: grp.add(h); new.append(h)
        frontier = new
    return sorted(grp)
def singer_matrix(poly=0b10011):
    cols = []
    for j in range(m):
        v = 1 << (j + 1)
        if v >> m: v ^= poly
        cols.append(v)
    return tuple(sum(((cols[j] >> i) & 1) << j for j in range(m)) for i in range(m))
def mpow(A, k):
    R = IDM
    for _ in range(k): R = matmul(A, R)
    return R
def make_H(spec):
    if spec == 'trivial': return [IDM]
    if spec.startswith('singer'):
        k = int(spec.split(':')[1]) if ':' in spec else 1
        return gen_group([mpow(singer_matrix(), k)])
    if spec.startswith('gens:'):
        gens = [tuple(int(v) for v in g.split('.')) for g in spec[5:].split('/')]
        return gen_group(gens)
    raise SystemExit('unknown H ' + spec)

class System:
    def __init__(s, H, W, dirs):
        s.H = H; n = s.n = len(H); s.W = W; s.M = M = len(dirs); s.dirs = dirs
        hi = {h: i for i, h in enumerate(H)}
        inv = [hi[next(x for x in H if matmul(x, h) == IDM)] for h in H]
        s.perm = np.array([[hi[matmul(H[g], H[k])] for k in range(n)] for g in range(n)])     # g.k
        s.iperm = np.array([[hi[matmul(H[inv[g]], H[k])] for k in range(n)] for g in range(n)])  # g^-1.k
        Pg = np.zeros((n, n, n))
        for g in range(n):
            for k in range(n): Pg[g, s.perm[g, k], k] = 1
        s.Pg = Pg
        # fibres: one Q per orbit of Q -> Q o h  (support -> h^-1 support)
        supp = {}
        for Q in range(NF): supp[frozenset(x for x in DIRS if QV[Q, x])] = Q
        seen = set(); reps = []
        for Q in range(1, NF):
            if Q in seen: continue
            S0 = frozenset(x for x in DIRS if QV[Q, x])
            orb = {supp[frozenset(apply(h, x) for x in S0)] for h in H}
            seen |= orb; reps.append(Q)
        s.reps = reps; F = s.F = len(reps)
        # Pi[f, t, k] = Q(k^-1 z_t);  O[f, k] = tr(Q o k^-1)
        kinv = [H[inv[k]] for k in range(n)]
        s.Pi = np.array([[[QV[Q, apply(kinv[k], z)] for k in range(n)] for z in dirs] for Q in reps], dtype=float)
        s.O = np.array([[sum(QV[Q, apply(kinv[k], 1 << i)] for i in range(m)) & 1 for k in range(n)] for Q in reps], dtype=float)
        s.tril = [(t, u) for t in range(M) for u in range(t)]
        s.nA = W * M * n; s.nB = M * W * n; s.nG = len(s.tril) * n; s.nu = s.nA + s.nB + s.nG
        s.OI = np.zeros((F, W * n, W * n))
        for l in range(W): s.OI[:, l * n:(l + 1) * n, l * n:(l + 1) * n] = np.einsum('fk,kj->fkj', s.O, np.eye(n))
    def unpack(s, x):
        n, W, M = s.n, s.W, s.M
        A = x[:s.nA].reshape(W, M, n); B = x[s.nA:s.nA + s.nB].reshape(M, W, n)
        G = np.zeros((M, M, n), dtype=complex)
        if s.nG:
            gv = x[s.nA + s.nB:].reshape(len(s.tril), n)
            ti = np.array(s.tril)
            G[ti[:, 0], ti[:, 1]] = gv
        return A, B, G
    def blocks(s, C):
        """C[..., g] -> left-regular matrices [..., n, n]"""
        return np.einsum('...g,gjk->...jk', C, s.Pg)
    def eval(s, x, jac=True):
        n, W, M, F = s.n, s.W, s.M, s.F
        A, B, G = s.unpack(x)
        At = s.blocks(A).transpose(0, 2, 1, 3).reshape(W * n, M * n)
        Bt = s.blocks(B).transpose(0, 2, 1, 3).reshape(M * n, W * n)
        Gt = s.blocks(G).transpose(0, 2, 1, 3).reshape(M * n, M * n)
        pi = s.Pi.reshape(F, M * n)
        L = np.eye(M * n)[None] - pi[:, :, None] * Gt[None]
        Y = np.linalg.solve(L, pi[:, :, None] * Bt[None])                 # (F, Mn, Wn)
        R = np.einsum('ab,fbc->fac', At, Y) - s.OI
        if not jac: return R
        # Z = At L^-1  -> solve L^T Z^T = At^T
        Z = np.linalg.solve(L.transpose(0, 2, 1), np.broadcast_to(At.T, (F, M * n, W * n))).transpose(0, 2, 1)
        V = Z * pi[:, None, :]                                             # (F, Wn, Mn)
        Yb = Y.reshape(F, M, n, W * n)
        Vb = V.reshape(F, W * n, M, n).transpose(0, 2, 1, 3)              # (F, M, Wn, n)
        VbP = np.stack([Vb[..., s.perm[g]] for g in range(n)])            # (g, F, M, Wn, n): [:, k] = V_t[:, g k]
        PY = np.stack([Yb[:, :, s.iperm[g], :] for g in range(n)])        # (g, F, M, n, Wn): rows j = Y[g^-1 j]
        Wn = W * n
        J = np.zeros((F, Wn, Wn, s.nu), dtype=complex)
        # A[l, t, g]: rows of block l += PY[g][t]
        col = 0
        for l in range(W):
            for t in range(M):
                for g in range(n):
                    J[:, l * n:(l + 1) * n, :, col] = PY[g, :, t]; col += 1
        for t in range(M):
            for l in range(W):
                for g in range(n):
                    J[:, :, l * n:(l + 1) * n, col] = VbP[g, :, t]; col += 1
        if s.nG:
            ti = np.array(s.tril)
            # (g, F, P, Wn, n) @ (F, P, n, Wn)
            Gj = np.einsum('gfpak,fpkb->fpgab', VbP[:, :, ti[:, 0]], Yb[:, ti[:, 1]])
            J[..., col:] = Gj.reshape(F, len(s.tril) * n, Wn, Wn).transpose(0, 2, 3, 1)
        return R, J.reshape(F * Wn * Wn, s.nu)

def plain_dirs(W):
    return [1 << j for l in range(W) for j in range(m)]

def lm(sysm, x0, iters=400, tol=1e-26, maxtime=600, verbose=False, lam0=1e-2, stall=1e-4):
    x = x0.copy(); R, J = sysm.eval(x); r = R.ravel(); cost = float(np.vdot(r, r).real)
    lam = lam0; t0 = time.time(); hist = []
    for it in range(iters):
        JH = J.conj().T
        JHJ = JH @ J; g = JH @ r
        d = np.real(np.diag(JHJ)).copy(); d[d < 1e-12] = 1e-12
        while True:
            try:
                dx = np.linalg.solve(JHJ + lam * np.diag(d), -g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = x + dx
            Rn = sysm.eval(xn, jac=False); rn = Rn.ravel(); cn = float(np.vdot(rn, rn).real)
            if cn < cost:
                lam = max(lam / 3, 1e-15); break
            lam *= 4
            if lam > 1e12: break
        if lam > 1e12: break
        x = xn; cost = cn
        R, J = sysm.eval(x); r = R.ravel()
        hist.append(cost)
        if verbose and it % 10 == 0: print('  it %d cost %.3e lam %.1e t %.0fs' % (it, cost, lam, time.time() - t0), flush=True)
        if cost < tol or time.time() - t0 > maxtime: break
        if it > 60 and hist[-50] - cost < stall * hist[-50] and cost > 1e-8: break
    return x, cost, it

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--H', default='trivial'); ap.add_argument('--W', type=int, default=1); ap.add_argument('--M', type=int, default=0)
    ap.add_argument('--dirs', default=''); ap.add_argument('--plain', action='store_true')
    ap.add_argument('--restarts', type=int, default=2); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--iters', type=int, default=400); ap.add_argument('--maxtime', type=float, default=600)
    ap.add_argument('--scale', type=float, default=0.5); ap.add_argument('--check', action='store_true')
    ap.add_argument('--save', default=''); ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    H = make_H(a.H)
    if a.plain: dirs = plain_dirs(a.W)
    elif a.dirs: dirs = [int(v) for v in a.dirs.split(',')]
    else: raise SystemExit('need --dirs or --plain')
    S = System(H, a.W, dirs)
    print('H=%s |H|=%d W=%d M=%d dirs=%s fibres=%d unknowns=%d' % (a.H, S.n, a.W, S.M, dirs, S.F, S.nu), flush=True)
    rng = np.random.default_rng(a.seed)
    if a.check:
        x = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * 0.4
        R, J = S.eval(x); dx = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)); eps = 1e-6
        R2 = S.eval(x + eps * dx, jac=False)
        fd = (R2 - R).ravel() / eps; an = J @ dx
        print('FD check: rel err %.2e' % (np.linalg.norm(fd - an) / np.linalg.norm(an))); return
    best = None
    for rs in range(a.restarts):
        x0 = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * a.scale
        t0 = time.time()
        x, cost, it = lm(S, x0, iters=a.iters, maxtime=a.maxtime, verbose=a.verbose)
        print('restart %d: cost %.6e after %d iterations, %.0f s' % (rs, cost, it, time.time() - t0), flush=True)
        if best is None or cost < best[0]: best = (cost, x)
    print('BEST %.6e' % best[0])
    if a.save: np.save(a.save, best[1])

if __name__ == '__main__':
    main()
