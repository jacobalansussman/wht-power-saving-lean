"""flat_check.py -- is the residual of a skeleton CONSTANT in the content coefficients (structural obstruction) or a local
minimum?  Sample random complex alpha (and alpha with some coefficients 0, +-1, +-1/2) and report min/max of res^2."""
import sys
from symsolve import *
m, W, H, shape, menu, seed, nsk = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
for s in range(nsk):
    sk, cost = build_skeleton(m, W, H, shape, menu, seed * 100 + s)
    if sk is None or not sk.setup_final(6000):
        continue
    prep_fast(sk)
    N = sk.nalpha
    vals = []
    rng = random.Random(s)
    for k in range(60):
        al = np.array([complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(N)])
        if k % 3 == 1:
            al = np.array([rng.choice([0, 1, -1, 0.5, -0.5, 1j, -1j]) for _ in range(N)], dtype=complex)
        r, _ = vp_eval(sk, al, [], want_J=False)
        if r is not None:
            vals.append(float(np.sum(np.abs(r) ** 2)))
    print(f"skel {s} cost {cost} unknowns {N} dims {[(a, b) for (a, b, E) in sk.final]}: res^2 over 60 random alpha: min {min(vals):.6f} max {max(vals):.6f}", flush=True)
