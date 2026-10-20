"""one65.py -- w8-bpgl14: the all30-mode obstruction in O O F F O's own characters, reduced to ONE constituent tau = (6,5), n = 1.
Rows t in Omega (kernel slot + the 4 reachable slots), kernel block of rank 14 (= dim Z'), so rows + K' span Omega and
A = Phi(Omega) must lie in v (x) C^5 + P (P = Phi(K')) for ONE v in C^6.  In the quotient C^30 / P (dim q): the 4-dim A/P must lie
in the column space of M(v) = Q [v e_1^T .. v e_5^T] (q x 5): bilinear in (v, mu).  Alternating least squares over many random
starts; prints the smallest residual found (0 => the single-type condition can be met)."""
import sys
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
f = np.load(sys.argv[1]); chi = int(sys.argv[2]); mode = sys.argv[3]; ty = int(sys.argv[4]); starts = int(sys.argv[5])
DS = int(f['DS']); n = 30 * DS; d0 = int(f['d0']); nt = int(f['ntypes']); Us = [f[f'U{t}'] for t in range(nt)]
free = list(f['free_d4']) if mode == 'd4' else list(f['Fr_all30'])
zi = [d0 * DS + c for c in range(DS)]; fi = [d * DS + c for d in free for c in range(DS)]
idx = np.array(zi + fi)
def null(M, tol=1e-10):
    U, s, Vh = np.linalg.svd(M, full_matrices=True); k = int((s > tol * max(1, s[0])).sum()); return Vh[k:].conj().T
Om = np.zeros((len(idx), n), dtype=complex); Om[np.arange(len(idx)), idx] = 1
Kp = Om[len(zi):]
U = Us[ty]; m, _, d = U.shape
def Phi(X): return np.stack([X @ Uk for Uk in U], axis=2).reshape(X.shape[0], d * m)   # entry (s, k) at s*m + k
A = Phi(Om); P = Phi(Kp)
uP, sP, vP = np.linalg.svd(P, full_matrices=True); rp = int((sP > 1e-10 * sP[0]).sum())
Qc = vP[rp:].conj().T                      # (d m) x q : coordinates mod P  (x -> x Qc)
Aq = A @ Qc; ua, sa, va = np.linalg.svd(Aq, full_matrices=False); ra = int((sa > 1e-10 * sa[0]).sum())
Ab = va[:ra]                               # ra x q : basis of A/P
q = Qc.shape[1]
print(f'[one65] type {ty} (d {d}, m {m}); dim A {np.linalg.matrix_rank(A, tol=1e-9)}, dim P {rp}, quotient q = {q}, dim A/P = {ra}', flush=True)
E = np.zeros((d, m, d * m)); 
for s in range(d):
    for k in range(m): E[s, k, s * m + k] = 1
EQ = np.einsum('skx,xq->skq', E, Qc)       # image of the unit matrix e_s e_k^T mod P
rng = np.random.default_rng(1); best = 9
for st in range(starts):
    v = rng.standard_normal(d) + 1j * rng.standard_normal(d); v /= np.linalg.norm(v)
    for it in range(300):
        Mv = np.einsum('s,skq->kq', v, EQ)                 # m x q : rows = images of v e_k^T
        # mu: Ab ~ mu Mv  (ra x m)
        mu = np.linalg.lstsq(Mv.T, Ab.T, rcond=None)[0].T
        r1 = np.linalg.norm(mu @ Mv - Ab)
        # v: Ab ~ sum_s v_s (mu @ EQ[s])  -> linear in v
        Bs = np.stack([(mu @ EQ[s]).reshape(-1) for s in range(d)], axis=1)
        v = np.linalg.lstsq(Bs, Ab.reshape(-1), rcond=None)[0]
        nv = np.linalg.norm(v); v /= nv
    Mv = np.einsum('s,skq->kq', v, EQ); mu = np.linalg.lstsq(Mv.T, Ab.T, rcond=None)[0].T
    res = np.linalg.norm(mu @ Mv - Ab) / np.linalg.norm(Ab)
    sv = np.linalg.svd(np.vstack([Mv, Ab]), compute_uv=False)
    best = min(best, res)
    if st < 5 or res < 1e-6: print(f'  start {st}: relative residual {res:.3e}; singular values of [M(v); A/P] ' + ' '.join('%.1e' % x for x in sv), flush=True)
print(f'[one65] best relative residual over {starts} starts: {best:.3e}')
# generic rank of M(v) = the image of v (x) C^5 modulo P: the maximum over v (lower semicontinuity)
rks = []
for _ in range(20):
    v = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    Mv = np.einsum('s,skq->kq', v, EQ); s_ = np.linalg.svd(Mv, compute_uv=False); rks.append((int((s_ > 1e-9 * s_[0]).sum()), s_))
print(f'[one65] rank of M(v) = dim of (v (x) C^{m} + P)/P for 20 random v: {sorted(set(r for r, _ in rks))}; singular values (one v): ' +
      ' '.join('%.1e' % x for x in rks[0][1]))
print(f'[one65] => the single-type condition needs dim >= dim A/P = {ra}; max over v is the generic value above')
