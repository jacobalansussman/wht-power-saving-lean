"""verify_kappa.py -- w8-bpgl14: end-to-end check of a min-rank-completion solution with the ACTUAL slot-module matrices.
From the final subspaces N_tau (mrc.py output) build rho_tau(kappa) mapping the eigen-row space of rho_tau(e_a) onto N_tau,
R(kappa) = sum_(tau,k) U_k rho_tau(kappa) U_k^* (an element of the image of C[B'_m]: same block on every copy of a type), and
solve w R(e_a) R(kappa) = T for the top rows w, T = the kernel rows of H on Omega.  O O F F O: two kappas, twists p and p2,
N_tau split into N1 (dim <= n1) + N2 (dim <= n2); solve [w1 w2] [R(e1 kappa1); R(e2 kappa2)] = T.
Reports: E_a U_k = U_k rho(e_a) (decomposition consistent with the twist), the residual of the solve, and T's blocks:
rank of the d0 block, size of the 15 other distance-4 blocks.  argv: seed psi sigma npz solution.npz family a [psi1 psi2]"""
import sys, os, time
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decomp_lib import load_twists
t0 = time.time()
seed, psi, sname, dec, solf, fam, a = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7])
psi1, psi2 = (int(sys.argv[8]), int(sys.argv[9])) if len(sys.argv) > 9 else (1, 3)
L, D, sigJ, DS, p, p2 = load_twists(seed, psi, sname.split(','))
f = np.load(dec); g = np.load(solf)
nt = int(f['ntypes']); Us = [f[f'U{t}'] for t in range(nt)]; V = [g[f'V{t}'] for t in range(nt)]
n = 30 * DS; d0 = int(f['d0']); Sx = list(f['Sx']); free = list(f['free_d4'])
chars = L.chars
def Emat(pp, aa): return sum(chars[aa][j] * D.R(pp[j], sigJ, DS) for j in range(7)) / 7
if fam == 'p': twists = [(p, a)]
else: twists = [(p, (a - psi1) % 7), (p2, (a - psi2) % 7)]
Es = [Emat(pp, aa) for pp, aa in twists]
# consistency of the decomposition with the twists
dmax = 0.0
for E in Es:
    for U in Us:
        r0 = U[0].conj().T @ E @ U[0]
        for Uk in U: dmax = max(dmax, np.abs(E @ Uk - Uk @ r0).max())
print(f'[verify] {sname} {fam} a/chi {a}: max |E U_k - U_k rho(e)| = {dmax:.1e}  {time.time() - t0:.0f}s', flush=True)
# kappa per twist
Ms = []
alloc = [0] * nt
for (pp, aa), E in zip(twists, Es):
    Rk = np.zeros((n, n), dtype=complex)
    for t in range(nt):
        U = Us[t]; dt = U.shape[2]
        rho_e = U[0].conj().T @ E @ U[0]
        Uu, s, Vh = np.linalg.svd(rho_e); k = int((s > 1e-8).sum())
        A = Vh[:k]                                         # orthonormal rows: row space of rho(e)
        Vt = V[t]; take = Vt[alloc[t]:alloc[t] + k]; alloc[t] += take.shape[0]
        K = A[:take.shape[0]].conj().T @ take if take.shape[0] else np.zeros((dt, dt))
        for Uk in U: Rk += Uk @ K @ Uk.conj().T
    Ms.append(E @ Rk)
print(f'[verify] subspaces used per type (dim N_tau, allocated): ' + ' '.join(f'{V[t].shape[0]}/{alloc[t]}' for t in range(nt)), flush=True)
# T = kernel rows of H on Omega (from N)
from numpy.linalg import svd, eigh
cols = []
for U, Vt in zip(Us, V):
    dt = U.shape[2]
    if Vt.shape[0]:
        _, _, Vh = svd(Vt, full_matrices=True); Vp = Vh[Vt.shape[0]:]
    else: Vp = np.eye(dt)
    for Uk in U: cols.append(Uk @ Vp.conj().T)
G = np.concatenate(cols, axis=1)
zi = [d0 * DS + c for c in range(DS)]; fi = [d * DS + c for d in free for c in range(DS)]
idx = np.array(zi + fi); Go = G[idx]
ev, E_ = eigh(Go @ Go.conj().T)
ker = E_[:, ev < 1e-10]
T = np.zeros((ker.shape[1], n), dtype=complex); T[:, idx] = ker.conj().T     # rows t = v^* (H v = 0  =>  |t G| = 0)
# check T's rows really are in N(N) (|t G| small)
print(f'[verify] kernel rows {ker.shape[1]}; max |t G| = {np.abs(T @ G).max():.1e}; next eigenvalue {ev[ker.shape[1]]:.1e}', flush=True)
Mst = np.vstack(Ms)
w, *_ = np.linalg.lstsq(Mst.T, T.T, rcond=None); w = w.T                     # w Mst = T
res = np.linalg.norm(w @ Mst - T) / np.linalg.norm(T)
TT = w @ Mst
zb = TT[:, zi]; sx = np.concatenate([TT[:, d * DS:(d + 1) * DS] for d in Sx], axis=1)
sv = np.linalg.svd(zb, compute_uv=False)
print(f'[verify] solve w [R(e kappa)] = T: relative residual {res:.1e}; |w| {np.linalg.norm(w):.1e}', flush=True)
print(f'[verify] reconstructed rows: d0 block singular values ' + ' '.join('%.1e' % x for x in sv) +
      f'; max |entry| on the 15 other distance-4 slots {np.abs(sx).max():.1e}', flush=True)
print('done %.0fs' % (time.time() - t0))
