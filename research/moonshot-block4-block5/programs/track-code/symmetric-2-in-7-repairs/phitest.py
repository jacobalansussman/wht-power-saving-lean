"""phitest.py -- the FRAME-SUM relaxation (augmentation C[B] -> C, a ring map): K x_0 in sum_s C[B] u_s implies
e_0 (x) delta_Kframe in span_C{ Phi(u_s) },  Phi([g] x_l) = e_{frame(g)} (x) e_l.  For a skeleton, minimise over the
content coefficients alpha (LM on the relaxed residual) and compare with the full residual."""
import sys
from symsolve import *
m, W, H, shape, menu, seed, nsk = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
for s in range(nsk):
    sk, cost = build_skeleton(m, W, H, shape, menu, seed * 100 + s)
    if sk is None:
        continue
    G = sk.G
    # frame index per (generator, term)
    fidx = {}; cols = []
    for i, T in enumerate(sk.gens):
        for j, (l, g) in enumerate(T):
            key = (G.frm(g), l)
            if key not in fidx:
                fidx[key] = len(fidx)
    Kkey = (G.frm(G.KID), 0)
    if Kkey not in fidx:
        print(f"skel {s}: no generator reaches the K frame of x_0"); continue
    nrow = len(fidx)
    rowmap = [np.array([fidx[(G.frm(g), l)] for (l, g) in T]) for T in sk.gens]
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    def relres(al):
        cv = sk.coefs(al)
        A = np.zeros((nrow, len(sk.gens)), dtype=complex)
        for i in range(len(sk.gens)):
            np.add.at(A[:, i], rowmap[i], cv[i])
        t = np.zeros(nrow, dtype=complex); t[fidx[Kkey]] = 1
        beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
        return A @ beta - t
    best = np.inf
    rng = random.Random(s)
    for st in range(3):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        def full(x):
            al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
        x, c = lm(lambda x: relres(full(x)), x0, iters=40) if nf else (x0, float(np.sum(np.abs(relres(np.ones(N, dtype=complex))) ** 2)))
        best = min(best, c)
    print(f"skel {s}: cost {cost}: frame-sum relaxation, array 0: min residual^2 = {best:.4f} (full problem gives (W-1)/W = {(W-1)/W:.4f} per array)", flush=True)
