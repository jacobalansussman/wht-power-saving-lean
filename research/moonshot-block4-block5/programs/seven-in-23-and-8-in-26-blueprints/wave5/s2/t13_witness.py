"""t13_witness.py (session 2) -- an explicit WITNESS for the survivor's relaxation pass (independent of decide's rank logic):
rebuild t10's skeleton, take a generic vector (gamma, mu) of decide's null space, and check directly, character by
character (chi outside E), that  sum_k gamma_k v_chi(item k of F1) - mu_chi e_K  lies in the free span (inputs, O1, O2 and
their 30 items each), with every mu_chi != 0.  Prints the residuals and |mu_chi|.  Also saves gamma, mu."""
import sys, os, runpy
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.dirname(here))
g = runpy.run_path(os.path.join(here, 't10_build.py'), run_name='t10')
E, O1, O2, F1, psi = g['E'], g['O1'], g['O2'], g['F1'], g['psi']
import numpy as np
free = [E.inp(0), O1, O2]; shared = [F1]
chis = [k for k in range(E.W) if k != psi]
fmenus = [E.menu_cells(u) for u in free]; smenu = E.menu_cells(F1)
# rebuild decide's joint matrix to get its null space
blocks = []; data = []
for chi in chis:
    vecs = []
    for u, mu_ in zip(free, fmenus):
        vecs.append(E.vchi(u, chi))
        for mc in mu_: vecs.append(E.menu_vchi(mc, chi))
    svecs = [E.menu_vchi(mc, chi) for mc in smenu]
    cells = sorted(set(k for v in vecs + svecs for k in v) | {E.iK}); ci = {k: i for i, k in enumerate(cells)}
    A = np.zeros((len(cells), len(vecs)), dtype=complex)
    for j, v in enumerate(vecs):
        for k, x in v.items(): A[ci[k], j] += x
    Bm = np.zeros((len(cells), 30), dtype=complex)
    for j, v in enumerate(svecs):
        for k, x in v.items(): Bm[ci[k], j] += x
    e = np.zeros(len(cells), dtype=complex); e[ci[E.iK]] = 1
    U, s, _ = np.linalg.svd(A, full_matrices=False); r = int((s > 1e-8 * s[0]).sum()); U = U[:, :r]
    Q = lambda X: X - U @ (U.conj().T @ X)
    blocks.append((Q(Bm), Q(e))); data.append((A, Bm, e))
M = np.zeros((sum(b.shape[0] for b, _ in blocks), 30 + len(chis)), dtype=complex); r0 = 0
for i, (QB, qe) in enumerate(blocks):
    M[r0:r0 + QB.shape[0], :30] = QB; M[r0:r0 + QB.shape[0], 30 + i] = -qe; r0 += QB.shape[0]
_, s, Vh = np.linalg.svd(M); rank = int((s > 1e-8 * s[0]).sum()); null = Vh[rank:].conj().T
z = null @ (np.random.default_rng(7).standard_normal(null.shape[1]) + 1j * np.random.default_rng(8).standard_normal(null.shape[1]))
gamma, mu = z[:30], z[30:]
print('null space dim', null.shape[1], '; |mu_chi| =', np.round(np.abs(mu) / np.max(np.abs(mu)), 4))
for i, chi in enumerate(chis):
    A, Bm, e = data[i]
    target = Bm @ gamma - mu[i] * e                 # must lie in the column span of A (the free vectors)
    c, *_ = np.linalg.lstsq(A, target, rcond=None)
    print(f'  chi {chi}: |top F1 part - mu e_K - free combination| = {np.linalg.norm(A @ c - target):.2e}  (|mu| = {abs(mu[i]):.3f}, |F1 part| = {np.linalg.norm(Bm @ gamma):.3f})')
np.save(os.path.join(here, 'witness_gamma_mu.npy'), z)
