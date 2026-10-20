#!/usr/bin/env python3
"""scan.py (t2-block4): many LM trials of lmh.py's exact system for one address group H and one shape W in M.
Each trial: a random sequence of direction ORBIT representatives (directions matter only up to H, e3-rmbuild
module form), kept only if it passes exact deletion in every fibre (|D_Q| >= W n_odd(Q)); then LM restarts.
Writes one JSON line per trial to --out (flushed), and saves the x of any trial with cost < --hit.
Usage: scan.py --H spec --W 3 --M 11 --trials 50 --seed 1 --out out.jsonl [--bias plain] [--maxtime 3000]"""
import argparse, json, time, sys
import numpy as np
sys.path.insert(0, '.')
import lmh
np.seterr(all='ignore')
ap = argparse.ArgumentParser()
ap.add_argument('--H', required=True); ap.add_argument('--W', type=int, required=True); ap.add_argument('--M', type=int, required=True)
ap.add_argument('--trials', type=int, default=20); ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--restarts', type=int, default=1); ap.add_argument('--iters', type=int, default=300)
ap.add_argument('--out', required=True); ap.add_argument('--hit', type=float, default=1e-12)
ap.add_argument('--bias', default='none'); ap.add_argument('--maxtime', type=float, default=3000)
ap.add_argument('--dirs', default='', help='fixed direction list (then trials = restarts on it)')
ap.add_argument('--scale', type=float, default=0.5)
a = ap.parse_args()
H = lmh.make_H(a.H); n = len(H)
orbs = []; seen = set()
for x in lmh.DIRS:
    if x in seen: continue
    o = sorted({lmh.apply(h, x) for h in H}); seen |= set(o); orbs.append(o)
reps = [o[0] for o in orbs]
rng = np.random.default_rng(a.seed)
t_start = time.time()
f = open(a.out, 'a')
def deletion_ok(S):
    nD = S.Pi.sum(axis=(1, 2)); nodd = S.O.sum(axis=1)
    return bool(np.all(nD >= a.W * nodd - 1e-9))
plainorb = sorted({next(i for i, o in enumerate(orbs) if (1 << j) in o) for j in range(4)})
for tr in range(a.trials):
    if time.time() - t_start > a.maxtime: break
    if a.dirs: dirs = [int(v) for v in a.dirs.split(',')]
    else:
        for _ in range(10000):
            if a.bias == 'plain' and rng.random() < 0.7:
                dirs = [reps[plainorb[rng.integers(len(plainorb))]] if rng.random() < 0.8 else reps[rng.integers(len(reps))] for _ in range(a.M)]
            else:
                dirs = [reps[rng.integers(len(reps))] for _ in range(a.M)]
            S = lmh.System(H, a.W, dirs)
            if deletion_ok(S): break
        else:
            print('no deletion-passing sequence found'); break
    S = lmh.System(H, a.W, dirs)
    if not deletion_ok(S):
        print('trial %d dirs %s FAILS deletion' % (tr, dirs), flush=True); continue
    best = None
    for rs in range(a.restarts):
        x0 = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * a.scale
        t0 = time.time()
        x, cost, it = lmh.lm(S, x0, iters=a.iters, maxtime=900)
        if best is None or cost < best[0]: best = (cost, x, it, time.time() - t0)
    rec = dict(H=a.H, W=a.W, M=a.M, dirs=dirs, cost=best[0], iters=best[2], secs=round(best[3], 1), seed=a.seed, trial=tr)
    f.write(json.dumps(rec) + '\n'); f.flush()
    print(json.dumps(rec), flush=True)
    if best[0] < a.hit:
        np.save('hit_%s_W%dM%d_s%d_t%d.npy' % (a.H.replace(':', '').replace('/', '_'), a.W, a.M, a.seed, tr), best[1])
        print('HIT', flush=True)
