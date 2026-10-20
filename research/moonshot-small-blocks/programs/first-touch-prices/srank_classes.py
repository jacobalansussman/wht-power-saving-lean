# srank_classes.py -- EXHAUSTIVE (over the left classes of D b, as in ms-mixed c10/c24: 14 engine, 1680 widened) ranks and
# stable ranks of the hook single-term defect maps A_b = pi_{Y'} rho(D1 b)|_{E_T}, Y = N, T = R and K.
import numpy as np, collections, sys
from vn2 import A_of, rank, rho, OIDX, NIDX, nullsp, TGT
from blk3 import rand_engine, rand_widened, D
def canon(M):
    rows = []
    for r in M:
        nz = np.nonzero(np.abs(r) > 1e-9)[0]; r = r / r[nz[0]]
        rows.append(tuple(np.round(r.real, 6) + 0.0) + tuple(np.round(r.imag, 6) + 0.0))
    return tuple(sorted(rows))
if __name__ == "__main__":
    rng = np.random.default_rng(7)
    if __name__ == "__main__":
      for game, gen, NS in (('engine', rand_engine, 4000), ('widened', rand_widened, 60000)):
        reps = {}
        for t in range(NS):
            b = gen(rng); k = canon(D[1] @ b)
            if k not in reps: reps[k] = b
        tab = collections.Counter()
        for k, b in reps.items():
            row = []
            for t in ('R', 'K'):
                A = A_of(b, t); s = np.linalg.svd(A, compute_uv=False); r = int((s > 1e-8 * max(1, s[0])).sum())
                row.append((r, round(float((s ** 2).sum() / s[0] ** 2), 3) if r else 0.0))
            tab[tuple(row)] += 1
        print(f'{game}: {len(reps)} classes; ((rank R, srank R), (rank K, srank K)) -> count:')
        for k, v in sorted(tab.items()): print('   ', k, v)
