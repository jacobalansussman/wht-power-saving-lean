"""diag.py -- where does the residual live?  For a skeleton, solve, then split the residual of array 0 by (input, distance
of the row term from the kernel frame)."""
import sys
from symsolve import *
m, W, H, shape, menu, seed, nsk = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
for s in range(nsk):
    sk, cost = build_skeleton(m, W, H, shape, menu, seed * 100 + s)
    if sk is None or not sk.setup_final(6000):
        continue
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    best = (np.inf, None)
    for st in range(2):
        x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
        if nf == 0:
            best = (0, x0); break
        x, c = vp_lm(sk, N, free_idx, x0)
        if c < best[0]: best = (c, x)
    al = np.ones(N, dtype=complex)
    if nf: al[free_idx] = best[1][:nf] + 1j * best[1][nf:]
    G = sk.G
    # rebuild rows of array 0's closure to know the terms
    tok = []; gd = []
    for i, T in enumerate(sk.gens):
        d = {}
        for j, t in enumerate(T):
            d[t] = len(tok) + 1; tok.append((i, j))
        gd.append(d)
    rows, cols, ent = G.closure(gd, (0, G.KID), 6000)
    inv_rows = {r: t for t, r in rows.items()}
    r, betas = sk.residual(al, want_beta=True) if hasattr(sk, 'residual') else (None, None)
    fl = flat(sk, al)
    nr, nc, R_, C_, F_ = sk.fin[0]
    A = np.zeros((nr, nc), dtype=complex); np.add.at(A, (R_, C_), fl[F_])
    t = np.zeros(nr, dtype=complex); t[0] = 1
    beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
    res = A @ beta - t
    byk = {}
    for ri in range(nr):
        (l, gid) = inv_rows[ri]
        d = G.dist(gid)
        byk[(l, d)] = byk.get((l, d), 0) + abs(res[ri]) ** 2
    nK = {}
    for ri in range(nr):
        (l, gid) = inv_rows[ri]
        if G.dist(gid) == 0:
            nK[l] = nK.get(l, 0) + 1
    print(f"skel {s} cost {cost} res^2(all arrays)={best[0]:.3f}; array 0 rows {nr} cols {nc}; K-frame rows per input {nK}; residual by (input, dist to K): " + ", ".join(f"{k}:{v:.3f}" for k, v in sorted(byk.items()) if v > 1e-9), flush=True)
