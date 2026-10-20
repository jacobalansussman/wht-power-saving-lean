"""stest2.py -- w8-bpgl14: the STRUCTURED-top test of stest.py, fast and for all characters in one run.
Top rows w = E_d1 R_sigma(A), A in C[B]: block at slot d is in Alg t_d (Alg = sigma(C[B_d1]), t_d = sigma(beta(d1, b_d))).
For O O F F O in character chi: T = w1 R(e1) R(kappa1) + w2 R(e2) R(kappa2) (readout X = I WLOG for a full-rank kernel block).
For FIXED random w1, w2 (random elements of Wsp, A1 and A2 independent) T is linear in the kappa blocks G_(i,tau) (n_i x d_tau):
solve T|_(15 other d4 slots) = 0, T|_d0 = I.  The design D has |J| DS rows; full row rank => solvable for that w.
Rank by the eigenvalues of D D^*; solution g = D^* (D D^*)^-1 target; then an END-TO-END rebuild with the actual R matrices.
O O F O F (one twist): the same with one kappa (expected NOT onto: fewer unknowns than equations).
argv: seed psi sigma(comma) dec.npz family(p|pp2) chis(comma) [psi1 psi2] [rowsfree]"""
import sys, os, time, random
import numpy as np, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from decomp_lib import *
t0 = time.time()
seed, psi, sname, dec, fam = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
chis = [int(x) for x in sys.argv[6].split(',')]
psi1, psi2 = (int(sys.argv[7]), int(sys.argv[8])) if len(sys.argv) > 8 else (1, 3)
rowsfree = len(sys.argv) > 9 and sys.argv[9] == 'rowsfree'
L, D, sigJ, DS, p, p2 = load_twists(seed, psi, sname.split(','))
f = np.load(dec); nt = int(f['ntypes']); Us = [f[f'U{t}'] for t in range(nt)]
n = 30 * DS; d0 = int(f['d0']); Sx = list(f['Sx'])
assert int(f['DS']) == DS
J = np.array([d * DS + c for d in Sx for c in range(DS)] + [d0 * DS + c for c in range(DS)])
# Wsp generators, computing only slot d1 of each letter
gd1, gd1i = D.g[D.d1], minv(D.g[D.d1])
def act_d1(b):
    w = mm(gd1, b); dd = D.di[L.cell(w)]; return dd, mm(w, D.gi[dd])
letters = [b for (b, s) in L.letters]
H = []; tt = [None] * 30
for b in letters:
    dd, beta = act_d1(b)
    if dd == D.d1: H.append(b)
    if tt[dd] is None: tt[dd] = sigJ(beta)
rr = random.Random(5)
while len(H) < 600: H.append(mm(rr.choice(H), rr.choice(H)))
Alg = alg_basis([sigJ(act_d1(b)[1]) for b in H], rr)
print(f'[stest2] {sname} {fam} chis {chis}{" ROWS FREE" if rowsfree else ""}: dim Alg {len(Alg)} of {DS * DS}; '
      f'coset reps for all 30 slots: {all(x is not None for x in tt)}  {time.time() - t0:.0f}s', flush=True)
nrng = np.random.default_rng(seed + 123)
def rand_w():
    w = np.zeros((DS, n), dtype=complex)
    if rowsfree:
        return nrng.standard_normal((DS, n)) + 1j * nrng.standard_normal((DS, n))
    for d in range(30):
        Z = sum((nrng.standard_normal() + 1j * nrng.standard_normal()) * A for A in Alg)
        w[:, d * DS:(d + 1) * DS] = Z @ tt[d]
    return w
Rp = [D.R(x, sigJ, DS) for x in p]; Rp2 = [D.R(x, sigJ, DS) for x in p2]
def Emat(R, a): return sum(L.chars[a][j] * R[j] for j in range(7)) / 7
target = np.zeros((DS, len(J)), dtype=complex); target[:, -DS:] = np.eye(DS)
for chi in chis:
    tw = [(Rp, chi)] if fam == 'p' else [(Rp, (chi - psi1) % 7), (Rp2, (chi - psi2) % 7)]
    Es = [Emat(R, a) for R, a in tw]
    ws = [rand_w() for _ in Es]; us = [w @ E for w, E in zip(ws, Es)]
    Arows = []
    for E in Es:
        row = []
        for U in Us:
            rho = U[0].conj().T @ E @ U[0]; _, s, Vh = np.linalg.svd(rho); row.append(Vh[:int((s > 1e-8).sum())])
        Arows.append(row)
    cols = []; index = []
    for i, u in enumerate(us):
        for t, U in enumerate(Us):
            A = Arows[i][t]; nA = A.shape[0]
            if nA == 0: continue
            dt = U.shape[2]
            Ls = np.stack([u @ Uk @ A.conj().T for Uk in U])            # m x DS x nA
            Ws = np.stack([Uk.conj().T[:, J] for Uk in U])              # m x dt x |J|
            cols.append(np.einsum('kir,kcj->ijrc', Ls, Ws).reshape(DS * len(J), nA * dt)); index.append((i, t, nA, dt))
    Dm = np.concatenate(cols, axis=1)
    G = Dm @ Dm.conj().T
    ev = np.linalg.eigvalsh(G)
    rk = int((ev > 1e-12 * ev[-1]).sum())
    if rk == G.shape[0]:
        g = Dm.conj().T @ np.linalg.solve(G, target.reshape(-1))
    else:
        g, *_ = np.linalg.lstsq(Dm, target.reshape(-1), rcond=1e-12)
    res = np.linalg.norm(Dm @ g - target.reshape(-1)) / np.linalg.norm(target)
    # end to end
    pos = 0; Rk = [np.zeros((n, n), dtype=complex) for _ in us]
    for (i, t, nA, dt) in index:
        Gm = g[pos:pos + nA * dt].reshape(nA, dt); pos += nA * dt
        K = Arows[i][t].conj().T @ Gm
        for Uk in Us[t]: Rk[i] += Uk @ K @ Uk.conj().T
    T = sum(w @ E @ R for w, E, R in zip(ws, Es, Rk))
    T2 = T[:, J]
    print(f'[stest2] chi {chi} (a = {[a for _, a in tw]}): design {Dm.shape[0]} x {Dm.shape[1]}, rank {rk} (eig of D D^*: min {ev[0]:.1e}, '
          f'max {ev[-1]:.1e}); residual {res:.1e}; END-TO-END T = sum w_i R(e_i) R(kappa_i): other d4 slots max {np.abs(T2[:, :-DS]).max():.1e}, '
          f'd0 block - I max {np.abs(T2[:, -DS:] - np.eye(DS)).max():.1e}; |kappa| max {max(np.abs(R).max() for R in Rk):.1e}  {time.time() - t0:.0f}s', flush=True)
print('done %.0fs' % (time.time() - t0))
