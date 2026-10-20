# ex-dimension, 14 Oct 2026.  Dimension counting on unit varieties, module (normal) form, in a representation rho.
#
# Normal form (pf-algebraic Theorem 8 / Lemma 7; meta-border bsearch.py): W live arrays, n paid steps; step j has rank
# r_j and output u_j = rho(W_{r_j}) c_j with content c_j = sum_z a_{j,z} z over z in {inputs E_1..E_W, u_1..u_{j-1}},
# a_{j,z} in A = rho(C[B]) (span of the free letters).  Output t: O_t = sum_z f_{t,z} z over all sources; target
# O_t = rho(T_t) E_t with T_t = K (kernel) or W_k (rank-k target).  Optionally inputs restricted to a B-invariant
# subspace U (pf-algebraic Lemma 10: a further relaxation): E_l = inclusion of U in slot l.
#
# Program map Phi: params (all a's and f's, complex) -> outputs (O_1..O_W), each a d x (u W) matrix.
# Reported:  P = #params, E = #equations (= W * d * u * W), naive expected dim P - E,
#            r = generic Jacobian rank of Phi (dimension of the closure of the image), fibre dim P - r,
#            dimR = generic dim of R_n = A-span of all sources (r >= W dimR), full = (dimR == d u W),
#            o = dimension of the symmetry orbit of the target tuple under (post-mix by M_W(A)) x (pre-mix by M_W(A_U)).
# Rigorous consequences (see NOTE.md):  if K-tuple is reachable (even as a limit) then r >= o  (orbit test);
#   r == E (dominant) => no rank method in rho excludes the shape; full => the shape is exactly feasible in rho.
import sys, itertools, time, json
import numpy as np, warnings
warnings.filterwarnings('ignore')

I4 = np.array([1, 1j, -1, -1j])
def popc(x): return bin(x).count('1')
def perm_mat(n, p):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x, p[x]] = 1
    return M

def letters(m, game='engine'):
    # copied from ../border/bsearch.py (engine letters typed from proofs/ladder/ladder_check.py); 'widened' adds a
    # Toffoli-type permutation (m >= 3) and a one-point phase i (all Z4 phase patterns with the affine maps).
    n = 1 << m
    free = []
    for j in range(m): free.append(perm_mat(n, [x ^ (1 << j) for x in range(n)]))
    for j in range(1, m):
        free.append(perm_mat(n, [(x & ~((1 << 0) | (1 << j))) | (((x >> j) & 1) << 0) | (((x >> 0) & 1) << j)
                                  for x in range(n)]))
    if m >= 2: free.append(perm_mat(n, [x ^ ((x & 1) << 1) for x in range(n)]))
    if m >= 3: free.append(perm_mat(n, [((x << 1) & (n - 1)) | (x >> (m - 1)) for x in range(n)]))
    for z in range(1, n): free.append(np.diag([I4[popc(z & x) % 4] for x in range(n)]))
    free.append(1j * np.eye(n))
    if game == 'widened':
        free.append(np.diag([1j if x == n - 1 else 1 for x in range(n)]))          # one-point phase (CS at m=2)
        if m >= 3: free.append(perm_mat(n, [x ^ (4 if (x & 3) == 3 else 0) for x in range(n)]))  # Toffoli
    ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
    def dirm(z): return ca * np.eye(n) + cb * perm_mat(n, [x ^ z for x in range(n)])
    Wr = [np.eye(n, dtype=complex)]
    for j in range(m): Wr.append(dirm(1 << j) @ Wr[-1])
    return free, Wr          # Wr[r] = dir e_r ... dir e_1 ; Wr[m] = kernel K

def ext_power(g, k):
    n = g.shape[0]; subs = list(itertools.combinations(range(n), k))
    M = np.zeros((len(subs), len(subs)), complex)
    for a, S in enumerate(subs):
        for b, T in enumerate(subs): M[a, b] = np.linalg.det(g[np.ix_(S, T)])
    return M

def rep_fn(name):
    def f(g):
        if name == 'V': return g
        if name == 'VV': return np.kron(g, g)
        if name == 'VVb': return np.kron(g, g.conj())
        if name == 'VVV': return np.kron(np.kron(g, g), g)
        if name.startswith('VL'): return np.kron(g, ext_power(g, int(name[2:])))
        if name.startswith('L'): return ext_power(g, int(name[1:]))
        raise ValueError(name)
    return f

def algebra_basis(gens, tol=1e-9, chunk=256):
    d = gens[0].shape[0]
    B = (np.eye(d, dtype=complex) / np.sqrt(d)).reshape(1, -1); frontier = B.copy()
    while frontier.shape[0]:
        newrows = []
        for g in gens:
            for c0 in range(0, frontier.shape[0], chunk):
                cand = (g @ frontier[c0:c0 + chunk].reshape(-1, d, d)).reshape(-1, d * d)
                for _ in range(2): cand = cand - (cand @ B.conj().T) @ B
                u, sv, vh = np.linalg.svd(cand, full_matrices=False)
                keep = sv > tol * max(1.0, np.sqrt(d))
                if keep.any():
                    nb = vh[keep]; nb = nb - (nb @ B.conj().T) @ B
                    qq, _ = np.linalg.qr(nb.T); nb = qq.T
                    B = np.concatenate([B, nb], 0); newrows.append(nb)
        frontier = np.concatenate(newrows, 0) if newrows else np.zeros((0, d * d), complex)
    return B.reshape(-1, d, d)

def module_closure(gens, X0, tol=1e-9):
    # dimension of the span of {w X0 : w a word in gens}, by closure in the space of d x u matrices (orthonormal basis)
    d, u = X0.shape
    B = X0.reshape(1, -1) / np.linalg.norm(X0); frontier = B.copy()
    while frontier.shape[0]:
        newrows = []
        for g in gens:
            cand = np.einsum('ij,ajk->aik', g, frontier.reshape(-1, d, u)).reshape(-1, d * u)
            for _ in range(2): cand = cand - (cand @ B.conj().T) @ B
            uu, sv, vh = np.linalg.svd(cand, full_matrices=False)
            keep = sv > tol
            if keep.any():
                nb = vh[keep]; nb = nb - (nb @ B.conj().T) @ B
                qq, _ = np.linalg.qr(nb.T); nb = qq.T
                B = np.concatenate([B, nb], 0); newrows.append(nb)
        frontier = np.concatenate(newrows, 0) if newrows else np.zeros((0, d * u), complex)
    return B.shape[0]

def btypes(F, d, rng, iters=400):
    # isotypic pieces of rho|B: eigenspaces of a random Hermitian element of the commutant (group average by
    # repeated generator averaging; converges because the identity is included).
    Y = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); Y = Y + Y.conj().T
    gens = F + [np.eye(d)]
    for _ in range(iters):
        Y = sum(g @ Y @ np.linalg.inv(g) for g in gens) / len(gens)
    Y = (Y + Y.conj().T) / 2
    w, V = np.linalg.eigh(Y)
    groups = []; cur = [0]
    for i in range(1, d):
        if abs(w[i] - w[cur[-1]]) < 1e-6 * max(1, abs(w).max()): cur.append(i)
        else: groups.append(cur); cur = [i]
    groups.append(cur)
    return [V[:, g] for g in groups], [float(w[g[0]]) for g in groups]

def rank_sv(M, rel=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    if s.size == 0 or s[0] == 0: return 0, None
    r = int((s > rel * s[0]).sum())
    gap = (float(s[r - 1] / s[0]), float(s[r] / s[0]) if r < s.size else 0.0)
    return r, gap

class Shape:
    def __init__(self, m, rep, game, W, steps, targets, U=None, seed=0):
        self.m = m; rng = np.random.default_rng(seed); self.rng = rng
        free, Wr = letters(m, game); f = rep_fn(rep)
        self.F = [f(g) for g in free]; self.Wr = [f(g) for g in Wr]
        self.d = d = self.F[0].shape[0]
        import os
        cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'Acache_m{m}_{rep}_{game}.npy')
        if os.path.exists(cache): self.A = np.load(cache)
        else:
            self.A = algebra_basis(self.F)
            if d >= 32: np.save(cache, self.A)
        self.s = self.A.shape[0]
        if U is None: self.U = np.eye(d, dtype=complex)
        elif U == 'diag':
            n0 = 1 << m; self.U = np.zeros((d, n0), complex)
            for x in range(n0): self.U[x * n0 + x, x] = 1
            for g in self.F: assert np.allclose(self.U @ (self.U.T @ (g @ self.U)), g @ self.U), 'U not B-invariant'
            self.pieces = None
        else:
            pieces, _ = btypes(self.F, d, rng)
            pieces.sort(key=lambda P: P.shape[1])
            self.U = pieces[U] if isinstance(U, int) else U
            self.pieces = [P.shape[1] for P in pieces]
        self.u = self.U.shape[1]; self.W = W; self.steps = steps; self.targets = targets
        self.gU = module_closure(self.F + [self.Wr[1]], self.U)          # dim C_G iota_U (span of all letter words on U)
        self.sG = None
        self.n = len(steps)
        self.slots = [W + j for j in range(self.n)] + [W + self.n] * W
        self.P = sum(self.slots) * self.s
        self.E = W * d * self.u * W
    def inputs(self):
        d, u, W = self.d, self.u, self.W
        Es = []
        for l in range(W):
            X = np.zeros((d, u * W), complex); X[:, l * u:(l + 1) * u] = self.U; Es.append(X)
        return Es
    def target(self, t):
        X = np.zeros((self.d, self.u * self.W), complex)
        X[:, t * self.u:(t + 1) * self.u] = self.Wr[self.targets[t]] @ self.U
        return X
    def random_point(self, scale=1.0):
        return [self.rng.normal(size=(k, self.s)) + 1j * self.rng.normal(size=(k, self.s)) for k in self.slots]
    def forward(self, coef, normalize=False):
        # normalize=True rescales each step output to norm 1 (same spans; avoids the dynamic range of deep programs)
        gens = self.inputs(); k = 0
        for j in range(self.n):
            M = np.tensordot(coef[k], self.A, axes=(1, 0)); k += 1
            c = sum(M[i] @ gens[i] for i in range(len(gens)))
            u = self.Wr[self.steps[j]] @ c
            gens.append(u / np.linalg.norm(u) if normalize else u)
        outs = []
        for t in range(self.W):
            M = np.tensordot(coef[k], self.A, axes=(1, 0)); k += 1
            outs.append(sum(M[i] @ gens[i] for i in range(len(gens))))
        return gens, outs
    def dimR(self, coef):
        gens, _ = self.forward(coef, normalize=True)
        vecs = np.concatenate([np.einsum('aij,gjk->agik', self.A, np.stack(gens)).reshape(-1, self.d * self.u * self.W)], 0)
        return rank_sv(vecs)[0]
    def jac_sketch(self, coef, k_cols, batch=None):
        if batch is None: batch = max(16, min(512, int(2e7 // (max(self.slots) * self.d * self.d))))
        # J @ S for a random Gaussian S (P x k_cols), by batched forward-mode derivatives (exact, holomorphic).
        d, uW = self.d, self.u * self.W
        gens0, _ = self.forward(coef)
        Ms = [np.tensordot(c, self.A, axes=(1, 0)) for c in coef]
        cols = []
        done = 0
        while done < k_cols:
            b = min(batch, k_cols - done); done += b
            def dMk(c):         # random tangent for one slot coefficient: (b, nsrc, d, d)
                S = (self.rng.normal(size=(b,) + c.shape) + 1j * self.rng.normal(size=(b,) + c.shape))
                return np.tensordot(S, self.A, axes=(2, 0))
            dg = [np.zeros((b, d, uW), complex) for _ in range(self.W)]
            kk = 0
            for j in range(self.n):
                M = Ms[kk]; D = dMk(coef[kk]); kk += 1
                dc = np.zeros((b, d, uW), complex)
                for i in range(len(dg)):
                    dc += D[:, i] @ gens0[i] + np.matmul(M[i], dg[i])
                del D
                dg.append(np.matmul(self.Wr[self.steps[j]], dc))
            outs = []
            for t in range(self.W):
                M = Ms[kk]; D = dMk(coef[kk]); kk += 1
                do = np.zeros((b, d, uW), complex)
                for i in range(len(dg)):
                    do += D[:, i] @ gens0[i] + np.matmul(M[i], dg[i])
                del D
                outs.append(do.reshape(b, -1))
            cols.append(np.concatenate(outs, 1).T)
        return np.concatenate(cols, 1)
    def orbit_dim(self):
        # tangent of {alpha T X}: entry (t,s) = alpha_ts T_s|U + T_t x_ts, alpha in A, x in A acting on U (A_U).
        AU = np.einsum('ij,ajk->aik', self.U.conj().T, np.einsum('aij,jk->aik', self.A, self.U))   # U* a U
        tot = 0
        for t in range(self.W):
            for s_ in range(self.W):
                Ts = self.Wr[self.targets[s_]] @ self.U; Tt = self.Wr[self.targets[t]] @ self.U
                v1 = np.einsum('aij,jk->aik', self.A, Ts).reshape(self.s, -1)
                v2 = np.einsum('ij,ajk->aik', Tt, AU).reshape(self.s, -1)
                tot += rank_sv(np.concatenate([v1, v2], 0))[0]
        return tot

def run(m, rep, game, W, steps, targets, U=None, seed=0, jac=True, label=''):
    t0 = time.time()
    S = Shape(m, rep, game, W, steps, targets, U, seed)
    coef = S.random_point()
    out = {'label': label, 'm': m, 'rep': rep, 'game': game, 'W': W, 'steps': steps, 'targets': targets,
           'cost': int(sum(steps)), 'd': S.d, 'u': S.u, 'dimA': S.s, 'P': S.P, 'E': S.E, 'naive': S.P - S.E,
           'dimF': S.d * S.u * W, 'dimCG': S.sG, 'Eeff': W * W * S.gU, 'naive_eff': S.P - W * W * S.gU}
    if U is not None: out['pieces'] = S.pieces
    out['dimR'] = S.dimR(coef); out['full'] = out['dimR'] == W * S.gU
    out['orbit'] = S.orbit_dim()
    if jac:
        k = min(S.P, W * W * S.gU) + 24
        J = S.jac_sketch(coef, k)
        r, gap = rank_sv(J)
        out['r'] = r; out['gap'] = gap; out['dominant'] = (r == out['Eeff']); out['fibre'] = S.P - r
        out['orbit_ok'] = out['orbit'] <= r
    out['sec'] = round(time.time() - t0, 1)
    return out

if __name__ == '__main__':
    spec = json.loads(sys.argv[1])
    print(json.dumps(run(**spec)), flush=True)
