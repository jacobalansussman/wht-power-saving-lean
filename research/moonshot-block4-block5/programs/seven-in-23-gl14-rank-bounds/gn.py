"""gn.py -- w8-bpgl14: Gauss-Newton / Levenberg-Marquardt for the minimum-rank completion (the alternating scheme of
mrc.py can stall on FEASIBLE instances: e.g. KEEP = {forced-zero type, (42,1,6)} at a = 0, r = 12, which is feasible by a
one-line argument, stalled at 3/16).
Rows t_i = base_i + D_i B (i = 1..r): B = basis of K' (= Omega' cap Free, Omega' = Omega cap kernels of the forced-zero types)
plus sections of the directions of Z' orthogonal to the current kernel block (a Grassmann chart; empty when r = dim Z').
Residual: for every type with n_tau < r m_tau, the normal part U_perp^* T_tau V_perp of the stacked coordinates T_tau
(r m_tau x d_tau) w.r.t. its rank-n_tau truncation; LM step on D; re-centre; repeat.  Residual 0 <=> feasible.
argv: npz family a mode r starts iters [psi1 psi2] [seed]"""
import sys, os, time
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
t0 = time.time()
f = np.load(sys.argv[1]); fam = sys.argv[2]; a = int(sys.argv[3]); mode = sys.argv[4]; r = int(sys.argv[5])
starts, iters = int(sys.argv[6]), int(sys.argv[7])
psi1, psi2 = (int(sys.argv[8]), int(sys.argv[9])) if len(sys.argv) > 9 else (1, 3)
rng = np.random.default_rng(int(sys.argv[10]) if len(sys.argv) > 10 else 1)
DS = int(f['DS']); n = 30 * DS; d0 = int(f['d0']); nt = int(f['ntypes'])
Us = [f[f'U{t}'] for t in range(nt)]
if fam == 'p': nn = [int(f[f'n{t}'][a]) for t in range(nt)]
else:
    a1, a2 = (a - psi1) % 7, (a - psi2) % 7
    nn = [min(Us[t].shape[2], int(f[f'n{t}'][a1]) + int(f[f'n2_{t}'][a2])) for t in range(nt)]
KEEP = os.environ.get('KEEP')
if KEEP is not None:
    keep = set(int(x) for x in KEEP.split(','))
    nn = [nn[t] if t in keep else Us[t].shape[2] for t in range(nt)]
free = list(f['free_d4']) if mode == 'd4' else list(f['Fr_all30'])
zi = [d0 * DS + c for c in range(DS)]; fi = [d * DS + c for d in free for c in range(DS)]
idx = np.array(zi + fi); nz = len(zi)
def null(M, tol=1e-10):
    U, s, Vh = np.linalg.svd(M, full_matrices=True); k = int((s > tol * max(1, s[0] if len(s) else 1)).sum())
    return Vh[k:].conj().T
zero = [t for t in range(nt) if nn[t] == 0]
if zero:
    Cm = np.concatenate([Us[t][k][idx] for t in zero for k in range(Us[t].shape[0])], axis=1)
    Bas = null(Cm.T)                         # columns: basis of Omega' (Omega coordinates); rows t = Bas[:, j]^T
else:
    Bas = np.eye(len(idx), dtype=complex)
Om = np.zeros((Bas.shape[1], n), dtype=complex); Om[:, idx] = Bas.T          # rows: basis of Omega' in C^n
Zp = Om[:, zi]                                                               # their d0 parts
u_, s_, vh_ = np.linalg.svd(Zp, full_matrices=False); kz = int((s_ > 1e-10).sum())
Kp = null(Zp.T).T @ Om                     # rows: basis of K' = Omega' cap Free  (combinations with zero d0 part)
q, _ = np.linalg.qr(Kp.T); Kp = q.T
# section of Z': rows of Omega' whose d0 parts are an orthonormal basis of Z'
coef = np.linalg.lstsq(Zp.T, vh_[:kz].T, rcond=None)[0].T                     # kz x dim Omega'
Sec = coef @ Om                                                              # kz rows, d0 parts = vh_[:kz]
Zbasis = vh_[:kz]
act = [t for t in range(nt) if nn[t] > 0 and nn[t] < r * Us[t].shape[0] and nn[t] < Us[t].shape[2]]
print(f'[gn] {sys.argv[1]} {fam} a/chi {a} mode {mode} r {r}: forced-zero types {zero}; dim Omega\' {Om.shape[0]}, dim Z\' {kz}, '
      f'dim K\' {Kp.shape[0]}; active types {act}; equations {sum((r * Us[t].shape[0] - nn[t]) * (Us[t].shape[2] - nn[t]) for t in act)}', flush=True)
if r > kz: print('[gn] r > dim Z\': impossible by the forced-zero bound'); sys.exit()
UK = {t: [Kp @ Uk for Uk in Us[t]] for t in act}          # coordinates of the K' basis in each copy: (dimK x d)
def coords(T, t): return np.concatenate([T @ Uk for Uk in Us[t]], axis=0)  # (m r) x d, row (k, i)
def resid(T):
    tot = 0.0; parts = []
    for t in act:
        C = coords(T, t); U, s, Vh = np.linalg.svd(C, full_matrices=True)
        k = nn[t]; Up, Vp = U[:, k:], Vh[k:].conj().T
        R = Up.conj().T @ C @ Vp
        parts.append((t, Up, Vp, R)); tot += np.linalg.norm(R) ** 2
    return tot, parts
best = None
for st in range(starts):
    # random start: Z_sp = random r-dim subspace of Z', random K' parts
    Q, _ = np.linalg.qr(rng.standard_normal((kz, r)) + 1j * rng.standard_normal((kz, r)))
    T = Q.T @ Sec + 0.3 * (rng.standard_normal((r, Kp.shape[0])) + 1j * rng.standard_normal((r, Kp.shape[0]))) @ Kp
    lam = 1e-2; obj, parts = resid(T)
    for it in range(iters):
        # variable directions: K' plus sections of Z' directions orthogonal to the current kernel block
        Zc = T[:, zi]; uq, sq, vq = np.linalg.svd(Zc, full_matrices=True)
        if kz > r:
            Cz = Zc @ Zbasis.conj().T                                               # r x kz (coordinates in Zbasis)
            Zperp = null(Cz).T.conj() if False else np.linalg.svd(Cz, full_matrices=True)[2][r:]   # complement in C^kz
        else:
            Zperp = np.zeros((0, kz))
        Bv = np.vstack([Kp, Zperp @ Sec]) if Zperp.shape[0] else Kp
        nb = Bv.shape[0]
        # Jacobian: d(U_perp^* C V_perp) / dD, C rows (k, i): dC[(k,i), :] = dD[i, :] (Bv U_k)
        Jb = []; rb = []
        for (t, Up, Vp, R) in parts:
            m = Us[t].shape[0]
            W = np.stack([Bv @ Uk @ Vp for Uk in Us[t]])            # m x nb x (d - n)
            Upr = Up.conj().T.reshape(-1, m, r)                       # (rm - n) x m x r   (row index of C = k r + i)
            J = np.einsum('pki,kjc->pcij', Upr, W).reshape(-1, r * nb)
            Jb.append(J); rb.append(R.reshape(-1))
        J = np.concatenate(Jb, axis=0); rv = np.concatenate(rb)
        JJ = J.conj().T @ J; g = J.conj().T @ rv
        while True:
            dD = -np.linalg.solve(JJ + lam * np.eye(JJ.shape[0]), g).reshape(r, nb)
            Tn = T + dD @ Bv
            # keep the d0 part orthonormal (left multiplication does not change feasibility)
            un, sn, vn = np.linalg.svd(Tn[:, zi], full_matrices=False)
            Tn = (un @ np.diag(1 / sn) @ un.conj().T) @ Tn
            on, pn = resid(Tn)
            if on < obj or lam > 1e8:
                break
            lam *= 10
        if on < obj: T, obj, parts = Tn, on, pn; lam = max(lam / 10, 1e-12)
        else: break
        if it % 10 == 0 or obj < 1e-26:
            print(f'    start {st} it {it}: objective {obj:.3e}  lam {lam:.0e}  {time.time() - t0:.0f}s', flush=True)
        if obj < 1e-26: break
    print(f'  start {st}: final objective {obj:.3e} after {it + 1} iterations; d0-block singular values ' +
          ' '.join('%.2f' % x for x in np.linalg.svd(T[:, zi], compute_uv=False)) + f'  {time.time() - t0:.0f}s', flush=True)
    if best is None or obj < best[0]: best = (obj, T)
obj, T = best
# report per-type numerical ranks of the final rows
line = []
for t in range(nt):
    C = coords(T, t); s = np.linalg.svd(C, compute_uv=False)
    line.append(f'({Us[t].shape[2]},{Us[t].shape[0]},{nn[t]}: s_n+1 {s[nn[t]] if nn[t] < len(s) else 0:.0e})')
print('[gn] BEST objective %.3e; per type ' % obj + ' '.join(line), flush=True)
np.savez_compressed(sys.argv[1].replace('.npz', '') + f'_gn_{fam}_{a}_{mode}_{r}.npz', T=T)
print('done %.0fs' % (time.time() - t0))
