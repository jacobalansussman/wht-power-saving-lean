"""mrc.py -- w8-bpgl14: the minimum-rank completion behind the exact kappa-free criterion (e3-bpexact NOTE s.4.4).
In one character, rows free, kappa free: the target T (r rows; Z on d0, 0 on the constrained slots, F free on the free
slots) is reachable iff for every B'_m-type tau the stacked coordinates [t_i U_k] (r m_tau x d_tau) have rank <= n_tau.
Equivalently: there are subspaces N_tau (dim n_tau) with every row t_i in N(N) = {t : t U_k in N_tau for all k, tau}.
For fixed N the best (Z, F) with Z Z^* = I comes from the Schur complement S(N) of H(N) = sum U_k (I - P_N) U_k^* on
Omega = (d0 | free): objective = sum of the r smallest eigenvalues of S(N).  Alternate with N_tau = top right singular
space of the stacked coordinates.  Objective 0 <=> a rank-r Z exists.
Bounds per type: 'p' (O O F O F: n_(tau,a)(p)), 'pp2' (O O F F O: n_(tau,a1)(p) + n_(tau,a2)(p2), a_i = chi - psi_i), or a
fixed list.  Multiple sigma blocks (joint decomposition): rows of block j live in block j's coordinates; r_j per block.
argv: npz family(p|pp2) a_or_chi mode(d4|all30) ranks(comma, one per sigma block) starts iters [psi1 psi2]"""
import sys, time
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)   # spurious FP flags of macOS Accelerate matmul
t0 = time.time()
EPS = 1e-10
_svd0 = np.linalg.svd
def _svd(a, full_matrices=True, compute_uv=True, **kw):
    # LAPACK gesdd occasionally fails to converge; retry on a tiny relative perturbation-free route (gesvd via scipy, or
    # a column-pivoted retry of gesdd on the conjugate transpose)
    try: return _svd0(a, full_matrices=full_matrices, compute_uv=compute_uv)
    except np.linalg.LinAlgError:
        try:
            import scipy.linalg as sl
            return sl.svd(a, full_matrices=full_matrices, compute_uv=compute_uv, lapack_driver='gesvd')
        except ImportError:
            r = _svd0(np.conj(a).T, full_matrices=full_matrices, compute_uv=compute_uv)
            if not compute_uv: return r
            U, s_, Vh = r; return np.conj(Vh).T, s_, np.conj(U).T
np.linalg.svd = _svd                      # Tikhonov term eps |F|^2: keeps F bounded; objective >= the true one
f = np.load(sys.argv[1]); fam = sys.argv[2]; a = int(sys.argv[3]); mode = sys.argv[4]
ranks = [int(x) for x in sys.argv[5].split(',')]; starts = int(sys.argv[6]); iters = int(sys.argv[7])
psi1, psi2 = (int(sys.argv[8]), int(sys.argv[9])) if len(sys.argv) > 9 else (1, 3)
DS = int(f['DS']); offs = f['offs']; nb = len(offs) - 1; n = 30 * DS; d0 = int(f['d0'])
nt = int(f['ntypes'])
Us = [f[f'U{t}'] for t in range(nt)]               # (m, n, d)
if fam == 'p':
    nn = [int(f[f'n{t}'][a]) for t in range(nt)]
elif fam == 'pp2':
    a1, a2 = (a - psi1) % 7, (a - psi2) % 7
    nn = [min(Us[t].shape[2], int(f[f'n{t}'][a1]) + int(f[f'n2_{t}'][a2])) for t in range(nt)]
else:
    raise ValueError(fam)
import os
KEEP = os.environ.get('KEEP')            # comma list of types to constrain; all others relaxed (n_tau = d_tau)
if KEEP is not None:
    keep = set(int(x) for x in KEEP.split(',') if x != '')
    nn = [nn[t] if t in keep else Us[t].shape[2] for t in range(nt)]
free = list(f['free_d4']) if mode == 'd4' else list(f['Fr_all30'])
cons = [d for d in range(30) if d != d0 and d not in free]
blocks = []
for j in range(nb):
    o, e = int(offs[j]), int(offs[j + 1])
    zi = np.array([d0 * DS + c for c in range(o, e)])
    fi = np.array([d * DS + c for d in free for c in range(o, e)])
    blocks.append((zi, fi))
print(f'[mrc] {sys.argv[1]} family {fam} a/chi {a} mode {mode} ranks {ranks}; types (d,m,n): ' +
      ' '.join(f'({U.shape[2]},{U.shape[0]},{k})' for U, k in zip(Us, nn)), flush=True)
rng = np.random.default_rng(int(time.time() * 1000) % 2 ** 31 if len(sys.argv) <= 10 else int(sys.argv[10]))

def Gmat(V):
    """H = G G^*, G = [U_k Vperp_tau^*] over all copies: t H t^* = |t G|^2 (stable form)"""
    cols = []
    for U, Vt in zip(Us, V):
        dt = U.shape[2]
        if Vt.shape[0]:
            _, _, Vh = np.linalg.svd(Vt, full_matrices=True); Vp = Vh[Vt.shape[0]:]
        else:
            Vp = np.eye(dt)
        for Uk in U: cols.append(Uk @ Vp.conj().T)
    return np.concatenate(cols, axis=1)

def Hmat(V):
    G = Gmat(V); return G @ G.conj().T

def zstep(G):
    """per block: minimise |t G|^2 over t in Omega with the d0-part orthonormal: z G_z (I - P_rowspace(G_f))"""
    rows = []; obj = 0.0; evs = []
    for (zi, fi), r in zip(blocks, ranks):
        Gz = G[zi]; Gf = G[fi]
        Uf, sf, Vfh = np.linalg.svd(Gf, full_matrices=False)
        k = int((sf > 1e-9 * max(1, sf[0])).sum())
        Vr = Vfh[:k]
        Gz2 = Gz - (Gz @ Vr.conj().T) @ Vr
        Sm = Gz2 @ Gz2.conj().T
        ev, E = np.linalg.eigh((Sm + Sm.conj().T) / 2)
        evs.append(ev)
        Z = E[:, :r].conj().T
        F = -((Z @ Gz @ Vr.conj().T) / sf[:k]) @ Uf[:, :k].conj().T
        T = np.zeros((r, n), dtype=complex); T[:, zi] = Z; T[:, fi] = F
        rows.append(T); obj += max(0.0, ev[:r].sum())
        obj_direct = np.linalg.norm(T @ G) ** 2
    return np.vstack(rows), obj, evs

def nstep(T):
    V = []
    for U, k in zip(Us, nn):
        C = np.concatenate([T @ Uk for Uk in U], axis=0)        # (m r) x d
        if k == 0: V.append(np.zeros((0, U.shape[2]))); continue
        _, s, Vh = np.linalg.svd(C, full_matrices=False)
        V.append(Vh[:k])
    return V

def check(V, T):
    """direct verification at the final N: dim(Omega cap N(N)) and its d0-projection rank per block; per-type ranks of T"""
    H = Hmat(V); out = []
    for (zi, fi), r in zip(blocks, ranks):
        idx = np.concatenate([zi, fi]); Ho = H[np.ix_(idx, idx)]
        ev, E = np.linalg.eigh(Ho)
        print('[check] smallest eigenvalues of H on Omega: ' + ' '.join('%.1e' % x for x in ev[:40]), flush=True)
        for thr in (1e-12, 1e-10, 1e-8):
            ker = E[:, ev < thr]
            zr = np.linalg.svd(ker[:len(zi)], compute_uv=False) if ker.shape[1] else np.zeros(1)
            out.append((thr, ker.shape[1], int((zr > 1e-4).sum()), ' '.join('%.0e' % x for x in zr[:16])))
    tr = []
    for U, k in zip(Us, nn):
        C = np.concatenate([T @ Uk for Uk in U], axis=0)
        s = np.linalg.svd(C, compute_uv=False)
        tr.append((k, float(s[k]) if k < len(s) else 0.0))
    return out, tr

best = None
for st in range(starts):
    V = []
    for U, k in zip(Us, nn):
        dt = U.shape[2]
        if k == 0: V.append(np.zeros((0, dt))); continue
        G = rng.standard_normal((dt, k)) + 1j * rng.standard_normal((dt, k)); q, _ = np.linalg.qr(G); V.append(q.T.conj()[:k].conj())
    hist = []
    for it in range(iters):
        T, obj, evs = zstep(Gmat(V))
        V = nstep(T)
        hist.append(obj)
        if obj < 1e-28: break
        if it > 300 and hist[-1] > 0.9999 * hist[-201]: break
        if it % 500 == 499: print(f'    it {it + 1}: objective {obj:.3e}', flush=True)
    if best is None or obj < best[0]: best = (obj, V, T, evs)
    print(f'  start {st}: iters {it + 1} objective {obj:.3e}  smallest S eigenvalues per block: ' +
          ' | '.join(' '.join('%.1e' % x for x in ev[:max(ranks) + 2]) for ev in evs) + f'  {time.time() - t0:.0f}s', flush=True)
obj, V, T, evs = best
out, tr = check(V, T)
print(f'[mrc] BEST objective {obj:.3e}; direct check per block (dim Omega cap N(N), d0-rank of it, next eigenvalue): {out}', flush=True)
print('[mrc] per type (n_tau, singular value number n_tau+1 of the stacked coordinates; ~0 = rank bound met): ' +
      ' '.join(f'({k},{s:.0e})' for k, s in tr), flush=True)
np.savez_compressed(sys.argv[1].replace('.npz', '') + f'_mrc_{fam}_{a}_{mode}_{"-".join(map(str, ranks))}.npz',
                    T=T, **{f'V{t}': v for t, v in enumerate(V)})
print('done %.0fs' % (time.time() - t0))
