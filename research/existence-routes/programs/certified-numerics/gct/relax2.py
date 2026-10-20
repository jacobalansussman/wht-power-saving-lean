# relax2.py -- W-in-M schemes in a QUOTIENT ("relaxation") of the block-3 engine game (a2-gct, 14 Oct 2026).
# Rewrite of ../relax.py after its gradient check failed (contents blew up along the gamma chains, the rank cut-off then
# dropped directions and the objective was not smooth). Changes: every paid-move output is NORMALISED (u_i -> u_i/|u_i|,
# which does not change the module), and the fixed rows coming from e_1..e_W are projected out exactly.
#
# Setting (proofs/ladder Lemma 0, module form). V = sum of irreducibles Sigma of Sp(6,2) (Pauli-blind quotient; matrices
# from meta/c2-finite-group/irreps.npz). R_Sigma = image of C[G] in End(V) (= full End of each irreducible, Burnside).
# Free algebra S_Sigma = image of C[P_0] = sum over P_0-types tau of End(tau), acting DIAGONALLY on all slots (copies)
# of tau in V. Paid letter h (one direction). Target k (the kernel). Node i: r_i = sum_l beta_il e_l +
# sum_{j<i} gamma_ij uhat_j, u_i = h r_i, uhat_i = u_i/|u_i| (beta, gamma in S_Sigma). Outputs are arbitrary
# S-combinations of e_l and u_i, so the scheme is right iff for every type tau the rows of the tau-view of k e_l lie in
# U_tau = rowspace(tau-views of e_l, u_1..u_M). Objective = sum_tau dist^2 / |target|^2 (variable projection).
# If the objective reaches 0 (or tends to 0 along a curve), T lies in the closure of the relaxed reachable set: then NO
# closed obstruction (rank, GCT occurrence/multiplicity, any polynomial) that only sees the irreducibles Sigma can
# exclude the W-in-M unit (relaxation duality, STATUS P4).
# Usage: relax2.py IRREPS(comma) W M MODE(target|plant) SEED [MAXITER] [INITSCALE]
import sys, json, time, os
import numpy as np
from scipy.optimize import minimize
NPZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'irreps_small.npz')
WORDS = None

def levi_words():
    g1 = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]); g2 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    seen = {np.eye(3, dtype=int).tobytes(): ()}; fr = [(np.eye(3, dtype=int), ())]
    while fr:
        nx = []
        for M, w in fr:
            for gi, g in enumerate((g1, g2)):
                N = (g @ M) % 2; kk = N.tobytes()
                if kk not in seen: seen[kk] = w + (gi,); nx.append((N, w + (gi,)))
        fr = nx
    return list(seen.values())

def load(js, rng):
    """direct sum of irreducibles js; returns n, h, k, types: list of (d_tau, [slot index arrays aligned])"""
    D = np.load(NPZ)
    global WORDS
    if WORDS is None: WORDS = sorted(levi_words(), key=len); assert len(WORDS) == 168
    blocks = []; off = 0; allpieces = []
    Hs, Ks = [], []
    for j in js:
        A, L, H = D[f'irr{j}_A'], D[f'irr{j}_L'], D[f'irr{j}_H']
        d = A.shape[1]
        Ael = []
        for bits in range(64):
            M = np.eye(d)
            for i in range(6):
                if bits >> i & 1: M = M @ A[i]
            Ael.append(M)
        Lel = {(): np.eye(d)}
        for w in WORDS:
            if w: Lel[w] = L[w[-1]] @ Lel[w[:-1]]
        Lel = list(Lel.values())
        S = rng.standard_normal((d, d)); S = S + S.T
        SL = sum(M @ S @ M.T for M in Lel) / 168; R = sum(M @ SL @ M.T for M in Ael) / 64
        w, V = np.linalg.eigh(R)
        cl = []; st = 0
        for i in range(1, d + 1):
            if i == d or w[i] - w[i - 1] > 1e-7 * max(1, abs(w).max()): cl.append((st, i)); st = i
        tests = [A[0], A[3], A[0] @ A[1], A[0] @ A[1] @ A[2], A[3] @ A[5], L[0], L[1], L[0] @ L[1], L[1] @ L[1] @ L[0],
                 A[0] @ L[1], A[3] @ L[0], A[2] @ A[4] @ L[1] @ L[0]]
        for (a, b) in cl:
            Q = V[:, a:b]
            sig = (b - a,) + tuple(round(float(np.trace(Q.T @ T @ Q)), 4) + 0.0 for T in tests)
            allpieces.append(dict(sig=sig, j=j, off=off, Q=Q, A=[Q.T @ x @ Q for x in A], L=[Q.T @ x @ Q for x in L]))
        Hs.append(H[0]); Ks.append(H[2]); blocks.append((off, d)); off += d
    n = off
    types = {}
    for p in allpieces: types.setdefault(p['sig'], []).append(p)
    Qg = np.zeros((n, n)); col = 0; tlist = []
    for sig in sorted(types, key=lambda s: (s[0], s)):
        ps = types[sig]; dt = sig[0]; ref = ps[0]; slots = []
        for p in ps:
            if p is ref: X = np.eye(dt)
            else:
                eqs = []
                for ga, gb in zip(ref['A'] + ref['L'], p['A'] + p['L']):
                    eqs.append(np.kron(np.eye(dt), ga.T) - np.kron(gb, np.eye(dt)))
                Mq = np.vstack(eqs); _, sv, vt = np.linalg.svd(Mq)
                assert sv[-1] < 1e-8 and (len(sv) < 2 or sv[-2] > 1e-6), ('intertwiner', sig, sv[-3:])
                X = vt[-1].reshape(dt, dt); X = X / np.sqrt(np.trace(X @ X.T) / dt)
                assert np.allclose(X @ X.T, np.eye(dt), atol=1e-8)
            Qp = np.zeros((n, dt)); Qp[p['off']:p['off'] + p['Q'].shape[0], :] = p['Q'] @ X
            Qg[:, col:col + dt] = Qp; slots.append(np.arange(col, col + dt)); col += dt
        tlist.append((dt, slots))
    assert col == n and np.allclose(Qg.T @ Qg, np.eye(n), atol=1e-8)
    Hb = np.zeros((n, n)); Kb = np.zeros((n, n))
    for (o, d), Hj, Kj in zip(blocks, Hs, Ks): Hb[o:o + d, o:o + d] = Hj; Kb[o:o + d, o:o + d] = Kj
    h = Qg.T @ Hb @ Qg; k = Qg.T @ Kb @ Qg
    assert np.allclose(h @ h.T, np.eye(n), atol=1e-8) and np.allclose(k @ k.T, np.eye(n), atol=1e-8)
    # the free algebra must commute with nothing in particular, but h and k must NOT be in it (sanity)
    return n, h, k, tlist


class Game:
    def __init__(self, n, h, k, tlist, W, M):
        self.n, self.h, self.k, self.T, self.W, self.M = n, h, k, tlist, W, M
        self.N = n * W
        self.sdim = sum(d * d for d, _ in tlist)
        self.npar = (M * W + M * (M - 1) // 2) * self.sdim
        E = np.zeros((W, n, self.N))
        for l in range(W): E[l, :, l * n:(l + 1) * n] = np.eye(n)
        self.E = E
        self.C = []   # orthonormal basis of the complement of rowspace(tau-views of e_1..e_W), per type
        for t in range(len(tlist)):
            Et = np.vstack([self.view(E[l], t) for l in range(W)])
            q, _ = np.linalg.qr(Et.T)                       # (mN) x (Wd), orthonormal columns
            P = np.eye(Et.shape[1]) - q @ q.T
            w, V = np.linalg.eigh(P)
            self.C.append(V[:, w > 0.5])
        self.slot_index = []
        for d, slots in tlist: self.slot_index.append(np.concatenate(slots))

    def view(self, x, t):            # tau-view: d_tau x (m_tau N)
        d, slots = self.T[t]
        return np.hstack([x[s, :] for s in slots])

    def unview(self, y, t, out):
        d, slots = self.T[t]; N = self.N
        for a, s in enumerate(slots): out[s, :] += y[:, a * N:(a + 1) * N]

    def smat(self, p):
        S = np.zeros((self.n, self.n)); o = 0
        for d, slots in self.T:
            blk = p[o:o + d * d].reshape(d, d); o += d * d
            for s in slots: S[np.ix_(s, s)] = blk
        return S

    def sgrad(self, G):
        out = []
        for d, slots in self.T:
            out.append(sum(G[np.ix_(s, s)] for s in slots).ravel())
        return np.concatenate(out)

    def unpack(self, x):
        sd = self.sdim; W, M = self.W, self.M; o = 0
        beta = [[None] * W for _ in range(M)]; gamma = {}
        for i in range(M):
            for l in range(W): beta[i][l] = self.smat(x[o:o + sd]); o += sd
        for i in range(M):
            for j in range(i): gamma[i, j] = self.smat(x[o:o + sd]); o += sd
        return beta, gamma

    def forward(self, x):
        beta, gamma = self.unpack(x); n, W, M = self.n, self.W, self.M
        Uh = []; nrm = []; Uraw = []
        for i in range(M):
            r = np.zeros((n, self.N))
            for l in range(W): r[:, l * n:(l + 1) * n] += beta[i][l]
            for j in range(i): r += gamma[i, j] @ Uh[j]
            u = self.h @ r; s = np.linalg.norm(u) + 1e-300
            Uraw.append(u); nrm.append(s); Uh.append(u / s)
        return beta, gamma, Uh, nrm

    def set_target(self, outs):
        self.Tt = [np.vstack([self.view(outs[l], t) for l in range(self.W)]) @ self.C[t] for t in range(len(self.T))]
        self.tnorm = sum(float((T * T).sum()) for T in self.Tt)

    def fg(self, x, need_grad=True, rtol=1e-11):
        beta, gamma, Uh, nrm = self.forward(x); n, W, M = self.n, self.W, self.M
        f = 0.0; dU = [np.zeros((n, self.N)) for _ in range(M)]; per = []
        for t, (d, slots) in enumerate(self.T):
            B = np.vstack([self.view(Uh[i], t) for i in range(M)]) @ self.C[t]
            T = self.Tt[t]
            if B.shape[1] == 0: per.append((0.0, 0, B.shape, 0.0)); continue
            u_, s, vt = np.linalg.svd(B, full_matrices=False)
            r = int((s > rtol * max(s[0], 1e-300)).sum())
            Vr = vt[:r]
            res = T - (T @ Vr.T) @ Vr; ft = float((res * res).sum()); f += ft
            per.append((ft, r, B.shape, float(s[r - 1] / s[0]) if r else 0.0))
            if need_grad and ft > 0:
                Bp = (vt[:r].T / s[:r]) @ u_[:, :r].T
                GB = (-2.0 * (T @ Bp).T @ res) @ self.C[t].T     # gradient w.r.t. the full tau-views of uhat
                for i in range(M):
                    self.unview(GB[i * d:(i + 1) * d, :], t, dU[i])
        if not need_grad: return f / self.tnorm, per
        gx = np.zeros(self.npar); sd = self.sdim
        Ubar = [g.copy() for g in dU]; dR = [None] * M
        for i in reversed(range(M)):
            # uhat_i = u_i / |u_i|: grad wrt u_i
            uh = Uh[i]; gu = (Ubar[i] - uh * float((uh * Ubar[i]).sum())) / nrm[i]
            dR[i] = self.h.T @ gu
            for j in range(i): Ubar[j] += gamma[i, j].T @ dR[i]
        o = 0
        for i in range(M):
            for l in range(W): gx[o:o + sd] = self.sgrad(dR[i][:, l * n:(l + 1) * n]); o += sd
        for i in range(M):
            for j in range(i): gx[o:o + sd] = self.sgrad(dR[i] @ Uh[j].T); o += sd
        return f / self.tnorm, gx / self.tnorm


def main():
    js = [int(a) for a in sys.argv[1].split(',')]; W = int(sys.argv[2]); M = int(sys.argv[3]); mode = sys.argv[4]
    seed = int(sys.argv[5]); maxit = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
    isc = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    rng = np.random.default_rng(seed); t0 = time.time()
    n, h, k, tl = load(js, np.random.default_rng(12345))
    g = Game(n, h, k, tl, W, M)
    tag = f'{"-".join(map(str, js))}_W{W}_M{M}_{mode}_s{seed}'
    print(f'[{tag}] irreps {js} n={n} types={[(d, len(s)) for d, s in tl]} dimS={g.sdim} W={W} M={M} '
          f'params={g.npar} ambient={W*W*n*n} complement dims={[c.shape[1] for c in g.C]} ({time.time()-t0:.1f}s)',
          flush=True)
    if mode == 'target':
        g.set_target([k @ g.E[l] for l in range(W)])
    else:   # plant: a random scheme's outputs (random S-combinations of e_l and u_i) -> must reach 0
        xs = rng.standard_normal(g.npar); _, _, Uh, _ = g.forward(xs)
        outs = []
        for l in range(W):
            o = sum(g.smat(rng.standard_normal(g.sdim)) @ g.E[ll] for ll in range(W))
            o = o + sum(g.smat(rng.standard_normal(g.sdim)) @ Uh[i] for i in range(M))
            outs.append(o)
        g.set_target(outs)
    x0 = rng.standard_normal(g.npar) * isc / np.sqrt(max(1, max(d for d, _ in tl)))
    f0, g0 = g.fg(x0); v = rng.standard_normal(g.npar); e = 1e-6
    fd = (g.fg(x0 + e * v, False)[0] - g.fg(x0 - e * v, False)[0]) / (2 * e)
    print(f'start f={f0:.6e} ranks={[(r, sh) for _, r, sh, _ in g.fg(x0, False)[1]]} gradcheck analytic={g0 @ v:.8e} '
          f'fd={fd:.8e} ({time.time()-t0:.1f}s)', flush=True)
    hist = []
    def cb(xk):
        cb.it += 1
        if cb.it % 100 == 0:
            fk, per = g.fg(xk, False)
            hist.append((cb.it, fk))
            print(f'  it {cb.it} f={fk:.4e} |x|={np.linalg.norm(xk):.3e} per-type(f,rank,minsv)='
                  f'{[(float("%.3e" % (a / g.tnorm)), r, float("%.1e" % c)) for a, r, _, c in per]} '
                  f'({time.time()-t0:.0f}s)', flush=True)
            np.save(f'sol_{tag}.npy', xk)
    cb.it = 0
    res = minimize(g.fg, x0, jac=True, method='L-BFGS-B', callback=cb,
                   options=dict(maxiter=maxit, maxcor=50, ftol=1e-30, gtol=1e-18, maxfun=maxit * 3))
    fk, per = g.fg(res.x, False)
    print(f'END [{tag}] f={fk:.6e} nit={res.nit} msg={res.message} |x|={np.linalg.norm(res.x):.3e} '
          f'per-type(f,rank,shape,minsv)={[(a / g.tnorm, r, s, c) for a, r, s, c in per]} ({time.time()-t0:.0f}s)',
          flush=True)
    np.save(f'sol_{tag}.npy', res.x)
    print(json.dumps(dict(js=js, W=W, M=M, mode=mode, seed=seed, f=fk, nit=int(res.nit), hist=hist)), flush=True)


if __name__ == '__main__':
    main()
