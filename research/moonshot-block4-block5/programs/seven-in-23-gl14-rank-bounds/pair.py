"""pair.py -- w8-bpgl14: the exact JOINT-COPY bound.  If the rows' kernel parts span an r-dim subspace of Z', then rows + K' +
(dim Z' - r) more vectors span Omega', so  delta_tau <= m_tau n_tau + (dim Z' - r)  for every type:  r <= dim Z' - max(delta - m n).
Phi: Omega' -> (C^d)^m, t -> (t U_1, ..., t U_m); A = Phi(Omega'), P = Phi(K'), K' = Omega' cap Free.  With Z_sp = Z' forced,
the rows need  A  inside  N^m + P  for one N of dim n.  Necessary: delta = dim A - dim P <= m n.  Prints the dims and the
structure of A/P (as a subspace of (C^d)^m / P).  argv: npz family a mode type"""
import sys
import numpy as np, warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
f = np.load(sys.argv[1]); fam, a, mode = sys.argv[2], int(sys.argv[3]), sys.argv[4]
DS = int(f['DS']); n = 30 * DS; d0 = int(f['d0']); nt = int(f['ntypes']); Us = [f[f'U{t}'] for t in range(nt)]
if fam == 'p': nn = [int(f[f'n{t}'][a]) for t in range(nt)]
else: nn = [min(Us[t].shape[2], int(f[f'n{t}'][(a - 1) % 7]) + int(f[f'n2_{t}'][(a - 3) % 7])) for t in range(nt)]
free = list(f['free_d4']) if mode == 'd4' else list(f['Fr_all30'])
zi = [d0 * DS + c for c in range(DS)]; fi = [d * DS + c for d in free for c in range(DS)]
idx = np.array(zi + fi)
def null(M, tol=1e-10):
    U, s, Vh = np.linalg.svd(M, full_matrices=True); k = int((s > tol * max(1, s[0])).sum()); return Vh[k:].conj().T
def rank(M, tol=1e-9):
    s = np.linalg.svd(M, compute_uv=False); return int((s > tol * max(1, s[0])).sum()), s
zero = [t for t in range(nt) if nn[t] == 0]
if zero:
    Cm = np.concatenate([Us[t][k][idx] for t in zero for k in range(Us[t].shape[0])], axis=1); Bas = null(Cm.T)
else:
    Bas = np.eye(len(idx), dtype=complex)
Om = np.zeros((Bas.shape[1], n), dtype=complex); Om[:, idx] = Bas.T
Kp = null(Om[:, zi].T).T @ Om
tys = range(nt) if sys.argv[5] == 'all' else [int(sys.argv[5])]
viol = []; excess = 0
for ty in tys:
    if nn[ty] == 0: continue
    U = Us[ty]; m, _, d = U.shape
    A = np.concatenate([Om @ Uk for Uk in U], axis=1)        # rows: Phi of the Omega' basis, in (C^d)^m
    P = np.concatenate([Kp @ Uk for Uk in U], axis=1)
    ra, sa = rank(A); rp, sp = rank(P)
    print(f'[pair] {fam} a/chi {a}: dim Z\' = {rank(Om[:, zi])[0]}; type {ty} (d {d}, m {m}, n {nn[ty]}), forced-zero types {zero}: dim Phi(Omega\') = {ra}, dim Phi(K\') = {rp}, '
          f'delta = {ra - rp}; need delta <= m n = {m * nn[ty]}')
    print('   singular values of Phi(K\') around the rank:', ' '.join('%.1e' % x for x in sp[max(0, rp - 3):rp + 3]))
    print('   singular values of Phi(Omega\') around the rank:', ' '.join('%.1e' % x for x in sa[max(0, ra - 3):ra + 3]))
    # structure of P inside (C^d)^m: is P = X^m for a subspace X of C^d (copy-symmetric)?  dims of P cap (C^d in copy k)
    for k in range(m):
        sel = np.zeros((m * d, d)); sel[k * d:(k + 1) * d] = np.eye(d)
        inter = rp + d - rank(np.vstack([P, sel.T]))[0]
        print(f'   dim P cap (copy {k} alone) = {inter}')

    if ra - rp > m * nn[ty]: viol.append(ty)
    excess = max(excess, ra - rp - m * nn[ty])
print(f'[pair] SUMMARY {fam} a/chi {a} {mode}: dim Z\' = {rank(Om[:, zi])[0]}; types violating delta <= m n (=> r = dim Z\' impossible): {viol}; JOINT-COPY BOUND r <= dim Z\' - max(delta - m n) = {rank(Om[:, zi])[0] - excess}')
