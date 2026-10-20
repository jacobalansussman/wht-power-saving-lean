#!/usr/bin/env python3
"""vscan.py (v-t2-block4): independent LM sampling (my solver vlm.py) for 3 in 11 at block 4, multi-target slice.
Mode 'perm': random orders of a fixed direction multiset; mode 'rand': random multisets of 11 directions that pass
deletion (rejection sampling), random order. Writes one JSON line per trial.
Usage: vscan.py mode targets multiset_or_- seed maxtime out.jsonl"""
import sys, json, time
import numpy as np
sys.path.insert(0, '.')
import vlm
mode, tgs, ms, seed, maxtime, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
tg = [vlm.parse(t) for t in tgs.split(',')]
rng = np.random.default_rng(seed); t0 = time.time(); k = 0
O = np.stack([vlm.pairing(vlm.Qs, S) for S in tg], 1).sum(1)
while time.time() - t0 < maxtime:
    if mode == 'perm':
        dirs = [int(v) for v in ms.split(',')]; dirs = [dirs[i] for i in rng.permutation(len(dirs))]
    else:
        for _ in range(200000):
            dirs = [int(v) for v in rng.integers(1, 16, 11)]
            if np.all(vlm.P[:, [z - 1 for z in dirs]].sum(1) >= O): break
        else:
            continue
    pr = vlm.Prob(tg, dirs)
    z0 = (rng.standard_normal(pr.n) + 1j * rng.standard_normal(pr.n)) * 0.5
    ts = time.time(); z, c, nf = vlm.solve(pr, z0, 2000)
    rec = dict(T=tgs, dirs=dirs, cost=c, norm=float(np.linalg.norm(z)), nfev=nf, secs=round(time.time() - ts, 1), seed=seed, k=k)
    open(out, 'a').write(json.dumps(rec) + '\n'); k += 1
    if c < 1e-12: np.save('vhit_%d_%d.npy' % (seed, k), z)
