# Special-point test, independently: get an exact (loss ~1e-29) 2-in-3 / 3-in-5 solution in V(x)V at block 2 with
# meta-border's optimiser (bsearch.Problem is used ONLY as an optimiser), then with my own float letters (typed
# independently, common.py conventions): (1) check every slot matrix lies in A = span rho(B); (2) re-run the program
# and check the outputs equal rho(K) E_t; (3) full (not sketched) Jacobian at the solution and at a random point, SVD.
import sys, os, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(1, os.path.join(HERE, '..', 'border'))
from scipy.optimize import minimize
W, T, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 0
def perm(n, f):
    M = np.zeros((n, n), complex)
    for x in range(n): M[x, f(x)] = 1
    return M
n = 4
free = [perm(n, lambda x: x ^ 1), perm(n, lambda x: x ^ 2), perm(n, lambda x: ((x & 1) << 1) | (x >> 1)),
        perm(n, lambda x: x ^ ((x & 1) << 1)), np.diag([1j if x & 1 else 1 for x in range(n)]),
        np.diag([1j if x & 2 else 1 for x in range(n)]), 1j * np.eye(n)]
ca, cb = (1 + 1j) / 2, (1 - 1j) / 2
D1 = ca * np.eye(n) + cb * perm(n, lambda x: x ^ 1); D2 = ca * np.eye(n) + cb * perm(n, lambda x: x ^ 2); K = D2 @ D1
key = lambda M: (np.round(M * 2, 6) + 0.0).tobytes()
seen = {key(np.eye(n, dtype=complex)): np.eye(n, dtype=complex)}; fr = [np.eye(n, dtype=complex)]
while fr:
    nf = []
    for M in fr:
        for g in free:
            X = g @ M; k = key(X)
            if k not in seen: seen[k] = X; nf.append(X)
    fr = nf
Bel = list(seen.values()); assert len(Bel) == 3072, len(Bel)
rho = lambda g: np.kron(g, g); d = 16
S = np.stack([rho(b).reshape(-1) for b in Bel]); u, sv, vh = np.linalg.svd(S, full_matrices=False)
rA = int((sv > 1e-9 * sv[0]).sum()); A = vh[:rA].conj().reshape(rA, d, d)   # orthonormal basis rows (conj so <A_a,X> = vdot)
Abasis = vh[:rA]                                                             # rows spanning vec(A)
print('dim A =', rA, flush=True)
rD, rK = rho(D1), rho(K)
import bsearch
Pb = bsearch.Problem(2, 'VV', W, T)
assert np.allclose(Pb.D, rD) and np.allclose(Pb.K, rK), 'bsearch letters differ from mine'
rng = np.random.default_rng(seed); best = None; t0 = time.time()
for start in range(8):
    x = rng.normal(size=Pb.npar) * 0.3
    res = minimize(Pb.loss_grad, x, jac=True, method='L-BFGS-B', options={'maxiter': 6000, 'maxfun': 24000, 'gtol': 1e-16, 'ftol': 1e-30})
    if best is None or res.fun < best.fun: best = res
    print('start', start, 'loss', res.fun, round(time.time() - t0), 's', flush=True)
    if res.fun < 1e-24: break
Mb = Pb.mats(Pb.coeffs(best.x))
# (1) membership in A (my basis)
mem = max(float(np.linalg.norm(M.reshape(-1) - (Abasis.conj().T @ (Abasis @ M.reshape(-1).conj())).conj()) / max(1e-300, np.linalg.norm(M))) for M in Mb)
# simpler membership: project vec(M) onto row space of Abasis
P_ = Abasis.T @ Abasis.conj()   # projector acting on column vectors vec(M): vec in rowspace(Abasis) => P_ v = v ? (rows are orthonormal)
mem = max(float(np.linalg.norm(Abasis.T @ (Abasis.conj() @ M.reshape(-1)) - M.reshape(-1)) / np.linalg.norm(M)) for M in Mb if np.linalg.norm(M) > 1e-12)
# (2) forward with my own code
slots = [W + j for j in range(T)] + [W + T] * W
Ms = []; k = 0
for ns in slots: Ms.append([Mb[k + i] for i in range(ns)]); k += ns
def inputs():
    out = []
    for l in range(W):
        X = np.zeros((d, d * W), complex); X[:, l * d:(l + 1) * d] = np.eye(d); out.append(X)
    return out
def forward(Ms):
    src = inputs(); k = 0
    for j in range(T):
        src.append(rD @ sum(Ms[k][i] @ src[i] for i in range(len(src)))); k += 1
    outs = []
    for t in range(W):
        outs.append(sum(Ms[k][i] @ src[i] for i in range(len(src)))); k += 1
    return src, outs
src, outs = forward(Ms)
resid = max(float(np.abs(outs[t] - rK @ inputs()[t]).max()) for t in range(W))
# (3) full Jacobian: one column per (slot, source, basis element a), directions = A_a (my orthonormal basis, unconjugated)
Adir = Abasis.reshape(rA, d, d)
def jac(Ms):
    src, _ = forward(Ms); P = sum(slots) * rA
    dsrc = [np.zeros((P, d, d * W), complex) for _ in range(W)]; k = 0; off = 0
    def block(kk, i, off):            # derivative contributions of slot kk, source i
        Z = np.zeros((P, d, d * W), complex); Z[off:off + rA] = np.matmul(Adir, src[i]); return Z
    for j in range(T):
        dc = np.zeros((P, d, d * W), complex)
        for i in range(len(dsrc)):
            dc[off:off + rA] += np.matmul(Adir, src[i]); off += rA
            dc += np.matmul(Ms[k][i], dsrc[i])
        k += 1; dsrc.append(np.matmul(rD, dc))
    cols = []
    for t in range(W):
        do = np.zeros((P, d, d * W), complex)
        for i in range(len(dsrc)):
            do[off:off + rA] += np.matmul(Adir, src[i]); off += rA
            do += np.matmul(Ms[k][i], dsrc[i])
        k += 1; cols.append(do.reshape(P, -1))
    J = np.concatenate(cols, 1); s = np.linalg.svd(J, compute_uv=False); s = s / s[0]
    r = int((s > 1e-9).sum()); return r, float(s[r - 1]), float(s[r]) if r < s.size else 0.0
r0 = jac(Ms)
Mr = [[(rng.normal(size=rA) + 1j * rng.normal(size=rA)) @ Abasis.reshape(rA, -1) for _ in range(ns)] for ns in slots]
Mr = [[v.reshape(d, d) for v in sl] for sl in Mr]
rg = jac(Mr)
print(json.dumps({'shape': f'{W} in {T}', 'loss': float(best.fun), 'max|slot matrix|': float(max(np.abs(M).max() for M in Mb)),
                  'membership_in_A_relerr': mem, 'residual_my_forward': resid, 'Eeff': W * W * 136,
                  'rank_at_solution(r, last kept sv, first dropped sv)': r0, 'generic_rank': rg}), flush=True)
