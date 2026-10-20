#!/usr/bin/env python3
"""clm.py (w8-rmexact, wave 8 resume): complex Levenberg-Marquardt with an ANALYTIC Jacobian for the exact no-phase
game C[Sym_m x| H], H = <sigma> cyclic of order n (Theorem A of e3-rmbuild; transfer form (**) of its Lemma C):

    for every binary quadratic form Q:   A Pi_Q (I + G Pi_Q)^-1 B  =  O_Q (x) I_W

A (W x M), B (M x W), G (M x M, strictly lower triangular in time) have entries in C[H], acting on the channel space
C^n left-regularly (entry (k', k) of the coefficient x is x[(k' - k) mod n]).  Pi_Q = diag over (move t, channel k) of
[Q(sigma^-k z_t) = 1]; O_Q = diag over channels of [sum_i Q(sigma^-k e_i) = 1] (target [I]).  One Q per H-orbit,
weighted by the orbit size.  Own code: shares nothing with tf.py / lm.py of e3-rmbuild.

The residual is holomorphic in the unknowns, so the Gauss-Newton matrix is J^H J (complex) and the step solves
(J^H J + mu diag(J^H J)) d = -J^H r.

Usage: clm.py --m 4 --gen singer:1 --W 2 --dirs 1,1,1,1,1,1,1 [--orbits 0] [--restarts 3] [--iters 200] ...
  --gen singer:p  = multiplication by alpha^p (alpha root of the standard primitive polynomial, polynomial basis)
  --gen c1,c2,..  = sigma(e_i) = c_i (integers, bit j = coordinate j)
  --dirs: the move directions in time order (integers), or 'u:a,b,..' = each a,b,.. repeated per --rep
"""
import argparse, json, time, sys, os
import numpy as np

PRIM = {2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101}


def lin_table(m, cols):
    N = 1 << m; t = [0] * N
    for v in range(N):
        r = 0
        for i in range(m):
            if (v >> i) & 1: r ^= cols[i]
        t[v] = r
    return t


def singer_cols(m, p):
    cols = []
    for i in range(m):
        v = 1 << i
        for _ in range(p):
            v <<= 1
            if (v >> m) & 1: v ^= PRIM[m]
        cols.append(v)
    return cols


def setup(m, cols, dirs):
    N = 1 << m
    sig = lin_table(m, cols)
    assert sorted(sig) == list(range(N)), "sigma not invertible"
    inv = [0] * N
    for v in range(N): inv[sig[v]] = v
    # order of sigma
    n = 1; cur = sig[:]
    while cur != list(range(N)):
        cur = [sig[v] for v in cur]; n += 1
    sk = [list(range(N))]                      # sk[k] = table of sigma^-k
    for k in range(1, n): sk.append([inv[v] for v in sk[-1]])
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    nf = 1 << len(pairs)
    def Qv(Q, x):
        s = 0
        for b, (i, j) in enumerate(pairs):
            if (Q >> b) & 1 and (x >> i) & 1 and (x >> j) & 1: s ^= 1
        return s
    val = np.array([[Qv(Q, x) for x in range(N)] for Q in range(nf)], dtype=np.int8)   # nf x N
    # orbits of forms under Q -> Q o sigma (as value tables)
    key = {val[Q].tobytes(): Q for Q in range(nf)}
    seen = np.zeros(nf, bool); reps = []
    for Q in range(nf):
        if seen[Q]: continue
        sz = 0; cur = Q
        while not seen[cur]:
            seen[cur] = True; sz += 1
            cur = key[val[cur][sig].tobytes()]
        reps.append((Q, sz))
    M = len(dirs); E = [1 << i for i in range(m)]
    PI = np.zeros((len(reps), M, n)); OD = np.zeros((len(reps), n)); WT = np.zeros(len(reps))
    for r, (Q, sz) in enumerate(reps):
        for k in range(n):
            for t, z in enumerate(dirs): PI[r, t, k] = val[Q, sk[k][z]]
            OD[r, k] = sum(int(val[Q, sk[k][e]]) for e in E) % 2
        WT[r] = sz
    return n, PI, OD, WT, reps


class Model:
    def __init__(s, W, M, n):
        s.W, s.M, s.n = W, M, n
        s.ti, s.tj = np.tril_indices(M, -1)
        s.nA = W * M * n; s.nB = M * W * n; s.nG = len(s.ti) * n
        s.P = s.nA + s.nB + s.nG
        a = np.arange(n); s.idx = (a[:, None] - a[None, :]) % n      # entry (k',k) -> coefficient index k'-k
    def unpack(s, th):
        W, M, n = s.W, s.M, s.n
        A = th[:s.nA].reshape(W, M, n); B = th[s.nA:s.nA + s.nB].reshape(M, W, n)
        G = np.zeros((M, M, n), complex); G[s.ti, s.tj] = th[s.nA + s.nB:].reshape(-1, n)
        return A, B, G
    def big(s, X):
        r, c, n = X.shape
        return X[:, :, s.idx].transpose(0, 2, 1, 3).reshape(r * n, c * n)
    def eval(s, th, PI, OD, WT, jac=False):
        """returns sum of weighted |r|^2, and (if jac) J^H J, J^H r accumulated over the given orbits."""
        W, M, n = s.W, s.M, s.n
        A, B, G = s.unpack(th)
        Ab, Bb, Gb = s.big(A), s.big(B), s.big(G)
        Imn = np.eye(M * n); Iw = np.eye(W)
        f = 0.0
        if jac:
            JHJ = np.zeros((s.P, s.P), complex); JHr = np.zeros(s.P, complex)
        for q in range(PI.shape[0]):
            p = PI[q].reshape(-1)                          # (t,k) flattened
            Nq = np.linalg.inv(Imn + Gb * p[None, :])      # (I + G Pi)^-1
            Y = Nq @ Bb                                     # Mn x Wn
            V = p[:, None] * Y                              # Pi N B
            U = (Ab * p[None, :]) @ Nq                      # A Pi N
            R = Ab @ V - np.kron(Iw, np.diag(OD[q]))
            w = WT[q]
            f += w * float(np.vdot(R, R).real)
            if not jac: continue
            sw = np.sqrt(w)
            # Jacobian columns, rows = flattened (row (l',k'), col (l2,k2)) of R
            Vt = V.reshape(M, n, W * n)                     # [t, k, c]
            Ut = U.reshape(W * n, M, n)                     # [r, t, k]
            # A-columns: (l,t,g): row-block l, rows k' = V[t, k'-g, :]
            sh = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n   # [k', g] -> k'-g
            VA = Vt[:, sh, :]                                # [t, k', g, c]
            JA = np.zeros((W, n, W * n, W, M, n), complex)
            for l in range(W):
                JA[l, :, :, l, :, :] = VA.transpose(1, 3, 0, 2)          # [k', c, t, g]
            JA = JA.reshape((W * n) * (W * n), s.nA)
            # B-columns: (t,l,g): col-block l, cols k = U[:, t, k+g]
            sh2 = (np.arange(n)[:, None] + np.arange(n)[None, :]) % n  # [k, g] -> k+g
            UB = Ut[:, :, sh2]                               # [r, t, k, g]
            JB = np.zeros((W * n, W, n, M, W, n), complex)
            for l in range(W):
                JB[:, l, :, :, l, :] = UB.transpose(0, 2, 1, 3)          # [r, k, t, g]
            JB = JB.reshape((W * n) * (W * n), s.nB)
            # G-columns: (t,t',g): - sum_k U[:, t, k+g] V[t', k, :]
            UBs = UB[:, s.ti]                                # [r, pair, k, g]
            Vs = Vt[s.tj]                                    # [pair, k, c]
            JG = -np.einsum('rpkg,pkc->rcpg', UBs, Vs, optimize=True).reshape((W * n) * (W * n), s.nG)
            J = np.concatenate([JA, JB, JG], axis=1) * sw
            r = R.reshape(-1) * sw
            JHJ += J.conj().T @ J
            JHr += J.conj().T @ r
        if jac: return f, JHJ, JHr
        return f
    def rvec(s, th, PI, OD, WT):
        W, M, n = s.W, s.M, s.n
        A, B, G = s.unpack(th); Ab, Bb, Gb = s.big(A), s.big(B), s.big(G)
        out = []
        for q in range(PI.shape[0]):
            p = PI[q].reshape(-1)
            Nq = np.linalg.inv(np.eye(M * n) + Gb * p[None, :])
            R = Ab @ (p[:, None] * (Nq @ Bb)) - np.kron(np.eye(W), np.diag(OD[q]))
            out.append(R.reshape(-1) * np.sqrt(WT[q]))
        return np.concatenate(out)


def plain_theta(mod, cols, m, W, noise, rng):
    """exact transfer data of the plain schedule (array l gets moves 4l..4l+3 in directions e_1..e_m, each a
    direction-1 move conjugated by sigma^k_j with sigma^k_j(e_1) = e_j), plus complex Gaussian noise."""
    n = mod.n; M = mod.M; assert M == m * W
    sig = lin_table(m, cols); kof = {}
    v = 1
    for k in range(n):
        kof.setdefault(v, k); v = sig[v]
    ks = [kof[1 << j] for j in range(m)] * W
    arr = [t // m for t in range(M)]
    A = np.zeros((W, M, n), complex); B = np.zeros((M, W, n), complex); G = np.zeros((M, M, n), complex)
    for t in range(M):
        A[arr[t], t, ks[t] % n] = 1; B[t, arr[t], (-ks[t]) % n] = 1
        for t2 in range(t):
            if arr[t2] == arr[t]: G[t, t2, (ks[t2] - ks[t]) % n] = 2
    th = np.concatenate([A.reshape(-1), B.reshape(-1), G[mod.ti, mod.tj].reshape(-1)])
    return th + noise * (rng.standard_normal(mod.P) + 1j * rng.standard_normal(mod.P)) / np.sqrt(2)


def pen_idx(mod, p):
    """indices in theta of A[:, p, :] and G[t, p, :] (t > p): the data by which move p reaches the outputs"""
    W, M, n = mod.W, mod.M, mod.n
    ia = [(l * M + p) * n + g for l in range(W) for g in range(n)]
    pairs = [i for i in range(len(mod.ti)) if mod.tj[i] == p]
    ig = [mod.nA + mod.nB + i * n + g for i in pairs for g in range(n)]
    return np.array(ia + ig)


def pen_idx_in(mod, p):
    """indices of B[p, :, :] and G[p, t, :] (t < p): the data feeding move p"""
    W, M, n = mod.W, mod.M, mod.n
    ib = [mod.nA + (p * W + l) * n + g for l in range(W) for g in range(n)]
    pairs = [i for i in range(len(mod.ti)) if mod.ti[i] == p]
    ig = [mod.nA + mod.nB + i * n + g for i in pairs for g in range(n)]
    return np.array(ib + ig, dtype=int)


class PenModel:
    """loss + lam |out|^2 + kappa |in|^2.  The per-move scalar gauge (out * c, in / c) leaves the program unchanged,
    so penalising only the outgoing data could be dodged by rescaling; with both terms the gauge-invariant product
    |out| |in| is what the penalty drives to zero (its minimum over the gauge is 2 sqrt(lam kappa) |out| |in|)."""
    def __init__(s, mod, idx, lam, idx_in=None, kappa=1.0):
        s.mod, s.idx, s.lam = mod, idx, lam; s.P = mod.P; s.n = mod.n
        s.idx_in = idx_in if idx_in is not None else np.array([], dtype=int); s.kappa = kappa
    def eval(s, th, PI, OD, WT, jac=False):
        out = s.mod.eval(th, PI, OD, WT, jac=jac)
        x = th[s.idx]; y = th[s.idx_in]
        pen = s.lam * float(np.vdot(x, x).real) + s.kappa * float(np.vdot(y, y).real)
        if not jac: return out + pen
        f, JHJ, JHr = out
        JHJ = JHJ.copy(); JHr = JHr.copy()
        JHJ[s.idx, s.idx] += s.lam; JHr[s.idx] += s.lam * x
        JHJ[s.idx_in, s.idx_in] += s.kappa; JHr[s.idx_in] += s.kappa * y
        return f + pen, JHJ, JHr


def lm(mod, th, PI, OD, WT, iters, mu=1e-2, tol=1e-26, verbose=False, maxtime=1e18):
    t0 = time.time()
    f, JHJ, JHr = mod.eval(th, PI, OD, WT, jac=True)
    for it in range(iters):
        dg = np.real(np.diag(JHJ)).copy() + 1e-12
        while True:
            Am = JHJ + mu * np.diag(dg)
            d = np.linalg.solve(Am, -JHr)
            th2 = th + d
            f2 = mod.eval(th2, PI, OD, WT)
            if f2 < f:
                th, mu = th2, max(mu / 3, 1e-15); break
            mu *= 4
            if mu > 1e14: return th, f, it
        f, JHJ, JHr = mod.eval(th, PI, OD, WT, jac=True)
        if verbose: print("   it %d f %.6e mu %.1e  (%.0fs)" % (it, f, mu, time.time() - t0), flush=True)
        if f < tol or time.time() - t0 > maxtime: return th, f, it
    return th, f, iters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, default=4); ap.add_argument('--gen', default='singer:1')
    ap.add_argument('--W', type=int, required=True); ap.add_argument('--dirs', required=True)
    ap.add_argument('--orbits', type=int, default=0); ap.add_argument('--restarts', type=int, default=3)
    ap.add_argument('--iters', type=int, default=200); ap.add_argument('--scale', type=float, default=0.5)
    ap.add_argument('--gscale', type=float, default=0.2); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--maxtime', type=float, default=1e18); ap.add_argument('--verbose', action='store_true')
    ap.add_argument('--label', default=''); ap.add_argument('--save', default='')
    ap.add_argument('--checkjac', action='store_true')
    ap.add_argument('--pattern', default='', help="init mask: for each move the arrays it touches, e.g. 0,0,0,1,1,1,01")
    ap.add_argument('--kill', type=int, default=-1, help='homotopy: penalise move p (its A column and outgoing Gamma) with --lams')
    ap.add_argument('--lams', default='1e-4,1e-3,1e-2,1e-1,1,10,100')
    ap.add_argument('--init', default='gauss', help='gauss | mono (each coefficient = random complex scalar times one random group element)')
    ap.add_argument('--stages', default='', help='orbit continuation: nested orbit subsets of these sizes (e.g. 4,8,16,72)')
    ap.add_argument('--shuffle', action='store_true', help='new random time order of the directions at every restart')
    ap.add_argument('--plain', type=float, default=-1.0, help='start from the exact plain schedule plus this much noise')
    a = ap.parse_args()
    m = a.m
    cols = singer_cols(m, int(a.gen.split(':')[1])) if a.gen.startswith('singer') else [int(v) for v in a.gen.split(',')]
    dirs = [int(v) for v in a.dirs.split(',')]
    n, PI, OD, WT, reps = setup(m, cols, dirs)
    M = len(dirs); mod = Model(a.W, M, n)
    rng0 = np.random.default_rng(a.seed + 999)
    if a.orbits and a.orbits < len(WT):
        keep = np.union1d(np.sort(rng0.choice(len(WT), size=a.orbits, replace=False)), [0])
    else:
        keep = np.arange(len(WT))
    print(json.dumps(dict(label=a.label, m=m, gen=cols, n=n, W=a.W, M=M, dirs=dirs, orbits_used=len(keep),
                          orbits_all=len(WT), unknowns_complex=mod.P)), flush=True)
    if a.checkjac:
        rng = np.random.default_rng(5); sl = keep[:3]
        th = (rng.standard_normal(mod.P) + 1j * rng.standard_normal(mod.P)) * 0.5
        d = rng.standard_normal(mod.P) + 1j * rng.standard_normal(mod.P)
        f, JHJ, JHr = mod.eval(th, PI[sl], OD[sl], WT[sl], jac=True)
        e = 1e-6
        rp = mod.rvec(th + e * d, PI[sl], OD[sl], WT[sl]); rm = mod.rvec(th - e * d, PI[sl], OD[sl], WT[sl])
        Jd = (rp - rm) / (2 * e)
        print('checkjac: |Jd|^2 fd %.10e  vs d^H JHJ d %.10e' % (np.vdot(Jd, Jd).real, np.vdot(d, JHJ @ d).real))
        r0 = mod.rvec(th, PI[sl], OD[sl], WT[sl])
        print('checkjac: <r, Jd> fd %s  vs (JHr)^H d %s' % (np.vdot(r0, Jd), np.vdot(JHr, d)))
        print('checkjac: f %.10e vs |r|^2 %.10e' % (f, np.vdot(r0, r0).real))
        return
    best = None
    for r in range(a.restarts):
        rng = np.random.default_rng(a.seed + r)
        if a.shuffle:
            dirs_r = [int(v) for v in rng.permutation(dirs)]
            n, PI, OD, WT, reps = setup(m, cols, dirs_r)
            print(json.dumps(dict(label=a.label, restart=r, dirs=dirs_r)), flush=True)
        th = (rng.standard_normal(mod.P) + 1j * rng.standard_normal(mod.P)) * a.scale / np.sqrt(2)
        th[mod.nA + mod.nB:] *= a.gscale / a.scale
        if a.init == 'mono':
            nco = mod.P // n
            T = np.zeros((nco, n), complex)
            T[np.arange(nco), rng.integers(0, n, nco)] = (rng.standard_normal(nco) + 1j * rng.standard_normal(nco)) / np.sqrt(2)
            th = T.reshape(-1) * a.scale
            th[mod.nA + mod.nB:] *= a.gscale / a.scale
        if a.pattern:
            touch = [set(int(c) for c in tok) for tok in a.pattern.split(',')]
            assert len(touch) == M
            Am, Bm, Gm = mod.unpack(th)
            for t in range(M):
                for l in range(a.W):
                    if l not in touch[t]: Am[l, t] = 0; Bm[t, l] = 0
                for t2 in range(t):
                    if not (touch[t] & touch[t2]): Gm[t, t2] = 0
            th = np.concatenate([Am.reshape(-1), Bm.reshape(-1), Gm[mod.ti, mod.tj].reshape(-1)])
        if a.plain >= 0:
            th = plain_theta(mod, cols, m, a.W, a.plain, rng)
            print(json.dumps(dict(label=a.label, restart=r, start='plain+noise', noise=a.plain,
                                  start_loss=mod.eval(th, PI[keep], OD[keep], WT[keep]))), flush=True)
        t0 = time.time()
        if a.kill >= 0:
            th, f, its = lm(mod, th, PI[keep], OD[keep], WT[keep], a.iters, maxtime=a.maxtime)
            idx = pen_idx(mod, a.kill)
            iin = pen_idx_in(mod, a.kill)
            print(json.dumps(dict(label=a.label, restart=r, phase='start', loss=f, movenorm=float(np.linalg.norm(th[idx])),
                                  innorm=float(np.linalg.norm(th[iin])))), flush=True)
            for lam in [float(v) for v in a.lams.split(',')]:
                pm = PenModel(mod, idx, lam, pen_idx_in(mod, a.kill), 1.0)
                th, fp, its = lm(pm, th, PI[keep], OD[keep], WT[keep], a.iters, tol=0.0, maxtime=a.maxtime)
                fm = mod.eval(th, PI[keep], OD[keep], WT[keep])
                print(json.dumps(dict(label=a.label, restart=r, lam=lam, loss=fm, movenorm=float(np.linalg.norm(th[idx])),
                                      innorm=float(np.linalg.norm(th[pen_idx_in(mod, a.kill)])),
                                      iters=its, secs=round(time.time() - t0, 1))), flush=True)
            th2 = th.copy(); th2[idx] = 0
            f = mod.eval(th2, PI[keep], OD[keep], WT[keep])
            print(json.dumps(dict(label=a.label, restart=r, phase='move removed', loss=f)), flush=True)
            th = th2
        elif a.stages:
            order = np.concatenate([[0], rng.permutation(np.arange(1, len(WT)))])
            for k in [int(v) for v in a.stages.split(',')]:
                sub = np.sort(order[:min(k, len(WT))])
                th, f, its = lm(mod, th, PI[sub], OD[sub], WT[sub], a.iters, verbose=a.verbose, maxtime=a.maxtime)
                print(json.dumps(dict(label=a.label, restart=r, stage=int(len(sub)), loss=f, iters=its,
                                      secs=round(time.time() - t0, 1))), flush=True)
                if f > 1e-12: break
            keep = sub
        else:
          th, f, its = lm(mod, th, PI[keep], OD[keep], WT[keep], a.iters, verbose=a.verbose, maxtime=a.maxtime)
        fall = mod.eval(th, PI, OD, WT)
        print(json.dumps(dict(label=a.label, restart=r, loss=f, loss_all=fall, iters=its,
                              secs=round(time.time() - t0, 1))), flush=True)
        if best is None or fall < best[0]: best = (fall, th)
    print("BEST", a.label, best[0], flush=True)
    if a.save: np.save(a.save, best[1])


if __name__ == '__main__':
    main()
