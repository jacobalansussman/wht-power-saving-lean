# v2_contphase.py -- v-ms-mixed: does ANY monomial letter with CONTINUOUS phases give a single-term entry D b that
# passes T2 for the two-row target R?  (The claim's scope line says "any monomial letters"; ms-mixed's c25 only sampled.)
# Classes up to left monomials: perfect matching {(x_k, y_k)} x torus q in U(1)^4 (phi_{x_k} = 1, phi_{y_k} = q_k).
# ||pi_O rho(D b(q)) E_R||_F^2 is a trigonometric polynomial in theta (q = e^{i theta}) with frequencies in {-2..2}^4;
# evaluated EXACTLY on a 32^4 grid by FFT, then the best grid points are refined by local minimisation.
import sys, os, itertools
import numpy as np
from scipy.optimize import minimize
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v1_common.py')).read())

GR = 32
results = []
for mi, m in enumerate(MAT):
    perm = [0] * 8
    for k, (x, y) in enumerate(m): perm[x] = 2 * k; perm[y] = 2 * k + 1
    sig, sgn, ex = rho_mono_sc(perm, [0] * 8)
    G0 = rhoD_c[:, sig] * sgn[None, :]          # rho(D) P_sigma S  (phases set to 1)
    # exponent vector alpha(j) for basis j = (x, T): alpha_k = [x == y_k] + [y_k in T]
    alpha = np.zeros((560, 4), int)
    for j, (x, T) in enumerate(BAS):
        for k, (xk, yk) in enumerate(m): alpha[j, k] = (x == yk) + (yk in T)
    # C_a = pi_O G0 Pi_a E_R   (E_R lives on the N coordinates)
    Cs = {}
    A0 = G0[np.ix_(OI, NI)]
    alN = alpha[NI]
    for a in set(map(tuple, alN)):
        msk = np.all(alN == np.array(a)[None, :], axis=1)
        Cs[a] = A0[:, msk] @ ER_c[msk, :]
    keys = list(Cs)
    Gd = np.zeros((GR,) * 4, complex)
    for a in keys:
        for b in keys:
            d = tuple((np.array(a) - np.array(b)) % GR)
            Gd[d] += np.vdot(Cs[b], Cs[a])      # tr(C_b^H C_a), coefficient of e^{i (a-b).theta}
    fgrid = np.real(np.fft.ifftn(Gd) * GR ** 4)   # f(theta_n) = sum_d G_d e^{+i d.theta_n}
    D_list = []; c_list = []
    for a in keys:
        for b in keys:
            D_list.append(np.array(a) - np.array(b)); c_list.append(np.vdot(Cs[b], Cs[a]))
    Dm = np.array(D_list); cv = np.array(c_list)
    def f(th):
        return float(np.real(np.sum(cv * np.exp(1j * (Dm @ th)))))
    def fdirect(th):
        q = np.exp(1j * th)
        A = sum(np.prod(q ** np.array(a)) * Cs[a] for a in keys)
        return float(np.real(np.vdot(A, A)))
    # check the trig-polynomial and FFT evaluations against the direct matrix evaluation
    for idx in ((3, 7, 11, 30), (0, 0, 0, 0), (8, 16, 24, 8)):
        th = 2 * np.pi * np.array(idx) / GR
        fd = fdirect(th)
        assert abs(fd - fgrid[idx]) < 1e-6 * max(1, fd) and abs(fd - f(th)) < 1e-6 * max(1, fd), (fd, fgrid[idx], f(th))
    th = np.random.default_rng(mi).uniform(0, 2 * np.pi, 4); fd = fdirect(th)
    assert abs(fd - f(th)) < 1e-6 * max(1, fd)
    gmin = fgrid.min()
    best = np.argsort(fgrid.ravel())[:8]
    loc = []
    for bi in best:
        th0 = 2 * np.pi * np.array(np.unravel_index(bi, fgrid.shape)) / GR
        r = minimize(f, th0, method='BFGS', options=dict(gtol=1e-12))
        loc.append(fdirect(r.x))
    results.append((min(loc), gmin, mi))
results.sort()
print("min over the torus of ||pi_O rho(D b(q)) E_R||_F^2, per matching (smallest 8):")
for v, g, mi in results[:8]: print(f"   matching {MAT[mi]}: local min {v:.6f} (grid min {g:.6f}, residual norm {np.sqrt(max(v,0)):.4f})")
print("overall smallest residual norm over all 105 matchings and all continuous phases (grid + local):",
      round(float(np.sqrt(max(results[0][0], 0))), 4))
