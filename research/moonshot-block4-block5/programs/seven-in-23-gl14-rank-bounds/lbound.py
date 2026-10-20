"""lbound.py -- w8-bpgl14: kappa-independent LOWER bounds for the minimum-rank completion (necessary conditions).
(1) types with n_tau = 0 force t U_k = 0 for every row (ubad.py's bound): Omega' = Omega cap (those kernels), Z' = its
    d0-projection, so r <= dim Z'.
(2) If r = dim Z', the row family must have d0-parts spanning Z', hence span(rows) + (Omega' cap Free) = Omega', and for
    every type: n_tau >= dim span{t U_k : t in Omega'} - dim span{t U_k : t in Omega' cap Free}.  A violation proves
    r <= dim Z' - 1.  (3) For r < dim Z' the single-type bound: with psi_(tau,k): Z' -> C^d / Q_tau,
    r <= dim Z' - rank(psi_(tau,1)) + n_tau  (any one copy k).
argv: npz family(p|pp2) a mode [psi1 psi2]"""
import sys
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
f = np.load(sys.argv[1]); fam = sys.argv[2]; a = int(sys.argv[3]); mode = sys.argv[4]
psi1, psi2 = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (1, 3)
DS = int(f['DS']); offs = f['offs']; nb = len(offs) - 1; n = 30 * DS; d0 = int(f['d0']); nt = int(f['ntypes'])
Us = [f[f'U{t}'] for t in range(nt)]
if fam == 'p': nn = [int(f[f'n{t}'][a]) for t in range(nt)]
else:
    a1, a2 = (a - psi1) % 7, (a - psi2) % 7
    nn = [min(Us[t].shape[2], int(f[f'n{t}'][a1]) + int(f[f'n2_{t}'][a2])) for t in range(nt)]
free = list(f['free_d4']) if mode == 'd4' else list(f['Fr_all30'])
def rank(M, tol=1e-9):
    if M.size == 0: return 0
    s = np.linalg.svd(M, compute_uv=False); return int((s > tol * max(1, s[0])).sum())
def null(M, tol=1e-9):
    U, s, Vh = np.linalg.svd(M, full_matrices=True); r = int((s > tol * max(1, s[0] if len(s) else 1)).sum())
    return Vh[r:].conj().T
for j in range(nb):
    o, e = int(offs[j]), int(offs[j + 1]); dsj = e - o
    zi = [d0 * DS + c for c in range(o, e)]; fi = [d * DS + c for d in free for c in range(o, e)]
    idx = np.array(zi + fi); nz = len(zi)
    # (1) forced-zero types
    zero = [t for t in range(nt) if nn[t] == 0]
    C0 = [Us[t][k][idx] for t in zero for k in range(Us[t].shape[0])]      # rows of Omega -> coordinates
    if C0:
        Cm = np.concatenate(C0, axis=1)                                      # |Omega| x (sum d m)
        Bas = null(Cm.T)                                                     # columns: basis of Omega' (in Omega coords)
    else:
        Bas = np.eye(len(idx))
    Zp = rank(Bas[:nz]); Kf = Bas @ null(Bas[:nz])                           # Omega' cap Free
    print(f'[lbound] block {j} ({dsj}): family {fam} a {a} mode {mode}; forced-zero types {[(Us[t].shape[2], Us[t].shape[0]) for t in zero]}; '
          f'dim Omega {len(idx)}, dim Omega\' {Bas.shape[1]}, dim Z\' = {Zp}, dim Omega\' cap Free {Kf.shape[1]}', flush=True)
    # (2) at r = dim Z'
    viol = []; line = []
    for t in range(nt):
        if nn[t] == 0: continue
        m, _, dt = Us[t].shape
        All = np.concatenate([(Bas.T @ Us[t][k][idx]) for k in range(m)], axis=0)    # coordinates of Omega' basis
        Fre = np.concatenate([(Kf.T @ Us[t][k][idx]) for k in range(m)], axis=0) if Kf.shape[1] else np.zeros((0, dt))
        ra, rf = rank(All), rank(Fre)
        lb = ra - rf
        line.append(f'({dt},{m},{nn[t]}: {ra}-{rf}={lb})')
        if lb > nn[t]: viol.append(t)
    print(f'[lbound]   at r = dim Z\' = {Zp}: per type (d,m,n: dim C(Omega\') - dim C(Omega\' cap Free) = bound): ' + ' '.join(line), flush=True)
    print(f'[lbound]   => r = {Zp} {"IMPOSSIBLE (types " + str(viol) + ")" if viol else "not excluded by this bound"}', flush=True)
    # (3) single-copy bound for r < dim Z'
    best = Zp
    for t in range(nt):
        if nn[t] == 0: continue
        m, _, dt = Us[t].shape
        Fre = np.concatenate([(Kf.T @ Us[t][k][idx]) for k in range(m)], axis=0) if Kf.shape[1] else np.zeros((0, dt))
        Q = null(Fre) if Fre.shape[0] else np.eye(dt)                  # columns: complement of rowspace(Fre)
        # psi_k: Omega' basis -> coordinates mod Q-span; restricted to a section of Z'
        for k in range(m):
            img = (Bas.T @ Us[t][k][idx]) @ Q                                  # (dim Omega') x (dt - rf)
            # kernel part along Free contributes nothing mod Q (by construction); rank over Omega' = rank over Z'
            rk = rank(img)
            bnd = Zp - rk + nn[t]
            best = min(best, bnd)
    print(f'[lbound]   single-copy bound: r <= {best}', flush=True)
