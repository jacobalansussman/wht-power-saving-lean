"""mt_lm.py (rm-theory): the MULTI-TARGET commutative slice, which by Theorem B is exactly the no-phase Singer game.

A multi-target slice program at block m: W arrays, array l owed [T_l] (T_l symmetric over F_2), moves [z z^T] in
time order z_1..z_M, scalar mixing free. Order-2 system (ex2-families Lemma 1 with one target per array; e3-slice6
Theorem A makes it exact): for every binary quadratic form Q with support D = {i : Q(z_i) = 1},
    A_D (I + H_D)^-1 B_D^T = diag_l(<Q, T_l>)          (H strictly upper triangular, refl.py's convention).
Trace placement targets: T_k = sigma^k sigma^kT (sigma = Singer cycle in self-dual coordinates); by Theorem B a
solution with these targets is a unit of C[Sym_m x| Z_n] with the same W and M.

Exact tests first (multi-target Lemma 2 / Theorem 1): deletion |D_Q| >= n_odd(Q); tight forms (equality, D nonempty);
tight-pair clique bound (pairwise joined live moves <= W). Then complex Levenberg-Marquardt (refl2.py's solver with the
target replaced by the diagonal of per-array parities). A residual reaching 1e-20 is EVIDENCE only; any hit must go to
the finite-field pipeline before it counts.

Usage: mt_lm.py m W_spec dirs_spec order restarts iters seed
  W_spec:  'k:0,1,4,...' (targets T_k), 'beta:1,2,3,...' (targets with diag = beta, field ints), 'I:4' (all I)
  dirs_spec: 'all' (every nonzero direction once), 'frames' (the union of the targets' frames), or 'list:...'
  order: 'lex', 'singer', or 'rand<seed>'
"""
import sys, time, itertools, random
import numpy as np
from gf import GF
from t_field import mons


def quad_forms(m):
    pairs = [(j, k) for j in range(m) for k in range(j, m)]
    nn = 1 << m
    zb = np.array([[(z >> j) & 1 for j in range(m)] for z in range(nn)], dtype=np.int64)
    mon = np.stack([zb[:, j] * zb[:, k] for (j, k) in pairs], axis=1)
    P = len(pairs)
    coeffs = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)
    vals = (coeffs @ mon.T) % 2
    return vals, coeffs, pairs


class Setup:
    def __init__(self, m, Wspec, dspec, order, seed=0):
        self.m = m
        F = self.F = GF(m); q, n = F.q, F.n
        fb = self.fb = F.self_dual_basis(1)
        co = {x: F.coords(x, fb) for x in range(q)}
        self.toint = {x: sum(c << j for j, c in enumerate(co[x])) for x in range(q)}   # field elt -> coord int
        def mat(a):
            cols = [co[F.mul(a, f)] for f in fb]
            return [[cols[j][i] for j in range(m)] for i in range(m)]
        def sym_of_frame(frame):
            S = [[0] * m for _ in range(m)]
            for g in frame:
                v = co[g]
                for i in range(m):
                    for j in range(m): S[i][j] ^= v[i] & v[j]
            return S
        # targets as (symmetric matrix, frame as field elements)
        kind, rest = Wspec.split(":")
        T = []
        if kind == "k":
            for k in map(int, rest.split(",")):
                frame = [F.mul(F.alpha(k), f) for f in fb]; T.append((sym_of_frame(frame), frame, f"T_k{k}"))
        elif kind == "beta":
            for b in map(int, rest.split(",")):
                # T with diagonal beta: frame gamma f_j with gamma^2 = beta
                g = F.pw(b, 1 << (m - 1))       # square root
                frame = [F.mul(g, f) for f in fb]; T.append((sym_of_frame(frame), frame, f"beta{b}"))
        elif kind == "gl":
            # explicit targets: indices into the list of all invertible non-alternating symmetric matrices
            allT = []
            ML_ = [(i, j) for i in range(m) for j in range(i, m)]
            from t_field import rank2
            for bits in range(1 << len(ML_)):
                S = [[0] * m for _ in range(m)]
                for idx, (i, j) in enumerate(ML_):
                    if (bits >> idx) & 1: S[i][j] = S[j][i] = 1
                if not any(S[i][i] for i in range(m)): continue
                if rank2([int("".join(map(str, r)), 2) for r in S], m) != m: continue
                allT.append(S)
            self.allT = allT
            for t in map(int, rest.split(",")):
                T.append((allT[t], None, f"gl{t}"))
        elif kind == "I":
            for _ in range(int(rest)):
                T.append((sym_of_frame(fb), list(fb), "I"))
        self.T = T; self.W = len(T)
        if dspec == "all":
            dirs = list(range(1, 1 << m))                     # coordinate ints
        elif dspec == "frames":
            dirs = [self.toint[g] for (_, fr, _) in T for g in fr]
        else:
            dirs = [int(x) for x in dspec.split(":")[1].split(",")]
        if order == "lex":
            dirs = sorted(dirs)
        elif order == "singer":
            # by discrete log of the field element with these coordinates
            inv = {v: x for x, v in self.toint.items()}
            dirs = sorted(dirs, key=lambda z: F.log[inv[z]])
        elif order.startswith("rand"):
            rng = random.Random(int(order[4:]) if order[4:] else seed); rng.shuffle(dirs)
        self.dirs = dirs; self.M = len(dirs)
        vals, coeffs, pairs = quad_forms(m)
        self.vals = vals[:, dirs]                              # F x M
        tvec = []
        for (S, _, _) in T:
            tvec.append(np.array([S[i][j] for (i, j) in pairs], dtype=np.int64))
        self.phi = (coeffs @ np.array(tvec).T) % 2             # F x W : parity of array l at form f

    def exact_tests(self):
        D = self.vals.sum(axis=1); nodd = self.phi.sum(axis=1)
        viol = int(np.sum(D < nodd))
        tight = np.nonzero((D == nodd) & (D > 0))[0]
        # vanishing odd forms: D = 0 but some array odd
        # clique: moves joined if they share a tight support; live = in some tight support
        M = self.M
        adj = np.zeros((M, M), dtype=bool); live = np.zeros(M, dtype=bool)
        for f in tight:
            S = np.nonzero(self.vals[f])[0]
            live[S] = True
            for a in S:
                adj[a, S] = True
        # max clique among live moves (M <= 31: simple branch and bound)
        nodes = [i for i in range(M) if live[i]]
        best = []
        def bk(R, P, X):
            nonlocal best
            if not P and not X:
                if len(R) > len(best): best = R[:]
                return
            if len(R) + len(P) <= len(best): return
            for v in list(P):
                bk(R + [v], [u for u in P if adj[v, u] and u != v], [u for u in X if adj[v, u] and u != v])
                P = [u for u in P if u != v]; X = X + [v]
        bk([], nodes, [])
        return dict(violations=viol, tight=len(tight), clique=len(best), W=self.W,
                    excluded=bool(viol > 0 or len(best) > self.W))


class MTSystem:
    """refl2.System2 with per-array targets."""
    def __init__(self, st):
        self.W, self.M = st.W, st.M
        W, M = self.W, self.M
        Ds, phis = [], []
        for f in range(st.vals.shape[0]):
            D = np.nonzero(st.vals[f])[0]
            if len(D) == 0:
                if st.phi[f].any(): raise ValueError("odd form vanishing on the moves")
                continue
            Ds.append(D); phis.append(st.phi[f].astype(float))
        self.Ds = Ds; self.phis = np.array(phis)
        self.nA = W * M
        iu = np.triu_indices(M, 1); self.iu = iu
        pidx = -np.ones((M, M), dtype=int); pidx[iu] = np.arange(len(iu[0]))
        self.nvar = 2 * self.nA + len(iu[0]); self.nres = len(Ds) * W * W
        self.groups = []
        sizes = np.array([len(D) for D in Ds])
        for d in sorted(set(sizes.tolist())):
            qq = np.nonzero(sizes == d)[0]
            Dg = np.array([Ds[i] for i in qq]); base = qq * W * W
            aa, bb = np.triu_indices(d, 1)
            hcols = 2 * self.nA + pidx[Dg[:, aa], Dg[:, bb]]
            self.groups.append((d, qq, Dg, base, aa, bb, hcols))

    def unpack(self, x):
        W, M = self.W, self.M
        A = x[:self.nA].reshape(W, M); B = x[self.nA:2 * self.nA].reshape(W, M)
        H = np.zeros((M, M), dtype=complex); H[self.iu] = x[2 * self.nA:]
        return A, B, H

    def residual(self, x, jac=False):
        W, M = self.W, self.M
        A, B, H = self.unpack(x)
        r = np.zeros(self.nres, dtype=complex)
        J = np.zeros((self.nres, self.nvar), dtype=complex) if jac else None
        pW = np.arange(W)
        for (d, qq, Dg, base, aa, bb, hcols) in self.groups:
            nf = len(qq)
            HD = H[Dg[:, :, None], Dg[:, None, :]]
            X = np.linalg.inv(np.eye(d)[None] + HD)
            AD = np.transpose(A[:, Dg], (1, 0, 2)); BD = np.transpose(B[:, Dg], (1, 0, 2))
            XB = X @ np.transpose(BD, (0, 2, 1)); AX = AD @ X
            Fm = AX @ np.transpose(BD, (0, 2, 1))
            R = Fm.copy()
            R[:, pW, pW] -= self.phis[qq]
            rows_all = base[:, None, None] + pW[None, :, None] * W + pW[None, None, :]
            r[rows_all.ravel()] = R.ravel()
            if jac:
                rows = base[:, None, None, None] + pW[None, :, None, None] * W + pW[None, None, None, :]
                rows = np.broadcast_to(rows, (nf, W, d, W))
                cols = np.broadcast_to(pW[None, :, None, None] * M + Dg[:, None, :, None], (nf, W, d, W))
                J[rows.ravel(), cols.ravel()] = np.broadcast_to(XB[:, None, :, :], (nf, W, d, W)).ravel()
                rows = base[:, None, None, None] + pW[None, :, None, None] * W + pW[None, None, :, None]
                rows = np.broadcast_to(rows, (nf, W, W, d))
                cols = np.broadcast_to(self.nA + pW[None, None, :, None] * M + Dg[:, None, None, :], (nf, W, W, d))
                J[rows.ravel(), cols.ravel()] = np.broadcast_to(AX[:, :, None, :], (nf, W, W, d)).ravel()
                if len(aa):
                    vals = -AX[:, :, None, aa] * np.transpose(XB[:, bb, :], (0, 2, 1))[:, None, :, :]
                    rows = np.broadcast_to(rows_all[:, :, :, None], vals.shape)
                    cols = np.broadcast_to(hcols[:, None, None, :], vals.shape)
                    J[rows.ravel(), cols.ravel()] = vals.ravel()
        return (r, J) if jac else r

    def solve(self, x0, iters=300, tol=1e-22, tmax=None):
        t0 = time.time(); x = x0.copy(); mu = 1e-3
        r, J = self.residual(x, jac=True); f = np.vdot(r, r).real
        for it in range(iters):
            JH = J.conj().T; g = JH @ r; JJ = JH @ J
            dg = np.real(np.diag(JJ)).copy()
            while True:
                step = np.linalg.solve(JJ + mu * np.diag(dg + 1e-12), -g)
                xn = x + step; rn = self.residual(xn); fn = np.vdot(rn, rn).real
                if fn < f:
                    x, f = xn, fn; mu = max(mu / 3, 1e-12); break
                mu *= 4
                if mu > 1e12: return x, f, it
            if f < tol or (tmax and time.time() - t0 > tmax): return x, f, it
            r, J = self.residual(x, jac=True)
        return x, f, iters


def main():
    m = int(sys.argv[1]); Wspec = sys.argv[2]; dspec = sys.argv[3]; order = sys.argv[4]
    restarts = int(sys.argv[5]); iters = int(sys.argv[6]); seed = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    tmax = float(sys.argv[8]) if len(sys.argv) > 8 else None
    st = Setup(m, Wspec, dspec, order, seed)
    et = st.exact_tests()
    print(f"m={m} targets={Wspec} dirs={dspec} order={order} W={st.W} M={st.M}: exact tests {et}", flush=True)
    if et["excluded"]:
        print("EXCLUDED by the exact tests; no search."); return
    sy = MTSystem(st)
    rng = np.random.default_rng(seed)
    best = None
    for rs in range(restarts):
        x0 = (rng.standard_normal(sy.nvar) + 1j * rng.standard_normal(sy.nvar)) / np.sqrt(2)
        x0[2 * sy.nA:] *= 0.3
        t0 = time.time()
        x, f, it = sy.solve(x0, iters=iters, tmax=tmax)
        print(f"  restart {rs}: residual {f:.3e} after {it} iterations, {time.time() - t0:.0f} s", flush=True)
        if best is None or f < best[1]: best = (x, f)
    if best is None: return
    np.save(f"mt_best_m{m}_{Wspec.replace(':', '_').replace(',', '-')}_{dspec[:6]}_{order}_s{seed}.npy", best[0])
    print(f"BEST residual {best[1]:.3e}")


if __name__ == "__main__":
    main()
