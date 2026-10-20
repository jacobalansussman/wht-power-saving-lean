# e3-boost enum.py: every small INTEGER relation on a register of m rows (m = 4..8) that the screen keeps:
# net vector v in [-K, K]^m (calls positive, tasks negative), primitive (gcd 1), with
#   frame count sum v >= 0, plain balance D(1) >= 0 (saves nothing alone), certificate A >= 0,
#   D(theta) < 0 somewhere in [0.9, 0.9983] (grid of 41 points),
#   then every rigorous certificate of certlib (hooks, lifts, projections, B, C_k, tensor scalings), exactly.
# Prints the survivors sorted by number of arrays (task weight), with their roots and tight certificates.
import sys, itertools
from math import gcd
from functools import reduce
import numpy as np
import certlib as C

def run(m, K, statuses=C.RIGOROUS + ("argued",), show=12):
    ds = np.arange(1, m + 1, dtype=float)
    ths = np.linspace(0.9, 0.9983, 41)
    PHI = ds[None, :] ** ths[:, None]                     # 41 x m
    Aw = np.array([2 ** m - 2 ** (m - d) for d in range(1, m + 1)], dtype=np.int64)
    vals = np.arange(-K, K + 1, dtype=np.int64)
    kept, n_window = [], 0
    # iterate over the first m-3 coordinates, vectorise the last 3
    tail = np.array(list(itertools.product(vals, repeat=3)), dtype=np.int64)
    for head in itertools.product(vals, repeat=m - 3):
        V = np.hstack([np.tile(np.array(head, dtype=np.int64), (len(tail), 1)), tail])
        ok = (V.sum(1) >= 0) & ((V * np.arange(1, m + 1)).sum(1) >= 0) & ((V * Aw).sum(1) >= 0)
        ok &= (V < 0).any(1) & (V > 0).any(1)
        if not ok.any(): continue
        Vk = V[ok]
        Dm = (Vk.astype(float)[:, None, :] * PHI[None, :, :]).sum(2).min(1)   # no BLAS (spurious FP warnings)
        Vk = Vk[Dm < -1e-12]
        n_window += len(Vk)
        for v in Vk:
            if reduce(gcd, [abs(int(x)) for x in v if x]) != 1: continue
            vd = {d + 1: int(x) for d, x in enumerate(v) if x}
            ver, bad, res = C.verdict(vd, m, statuses)
            if ver == "PASSES" or (ver == "excluded-weakly" and all(b[2] in ("sampled", "conj") for b in bad)):
                kept.append(vd)
    walk = [vd for vd in kept if C.walk_ok(vd, m)]
    kept.sort(key=lambda vd: (sum(-c for c in vd.values() if c < 0), sum(c for c in vd.values() if c > 0)))
    print(f"m={m} K={K}: window-negative primitive-or-not candidates after frame/plain/A prefilter: {n_window}; "
          f"primitive survivors of all rigorous certificates: {len(kept)}; of these, passing the walk condition: {len(walk)}")
    for vd in walk[:show]:
        tasks = {d: -c for d, c in vd.items() if c < 0}; calls = {d: c for d, c in vd.items() if c > 0}
        res = C.certify(vd, m, statuses)
        tight = sorted(set(n for (n, b, s, mg) in res if mg == 0 and b == 1))
        wmin = min(C.D(vd, 0.9 + (0.9983 - 0.9) * i / 400) for i in range(401))
        print(f"   tasks {tasks} calls {calls}: arrays {sum(tasks.values())}, D(1) = {C.D(vd, 1):.3g}, "
              f"min D on window {wmin:.3g}, roots {C.roots(vd)[:4]}, tight: {tight[:8]}")
    sys.stdout.flush()
    return kept, walk

if __name__ == "__main__":
    for a in sys.argv[1:]:
        m, K = map(int, a.split(":"))
        run(m, K, show=25)
