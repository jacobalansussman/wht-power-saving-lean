#!/usr/bin/env python3
"""scan_s.py (t2-block4): LM trials of the multi-target slice (lms.py) at block m for W arrays in M moves.
Targets: --targets fixed list, or --randT: T_1 = I and T_2..T_W random among the h h^T (all 420 at m = 4).
Directions: random sequences over all 2^m - 1 directions, kept only if they pass exact deletion in every character.
One JSON line per trial (flushed).  A trial with cost < --hit saves its x and is flagged HIT."""
import argparse, json, time, sys, itertools
import numpy as np
sys.path.insert(0, '.')
import lms
np.seterr(all='ignore')
ap = argparse.ArgumentParser()
ap.add_argument('--m', type=int, default=4); ap.add_argument('--W', type=int, required=True); ap.add_argument('--M', type=int, required=True)
ap.add_argument('--targets', default=''); ap.add_argument('--randT', action='store_true')
ap.add_argument('--trials', type=int, default=50); ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--iters', type=int, default=300); ap.add_argument('--out', required=True)
ap.add_argument('--hit', type=float, default=1e-12); ap.add_argument('--maxtime', type=float, default=3000)
ap.add_argument('--restarts', type=int, default=1); ap.add_argument('--dirs', default='')
ap.add_argument('--share', action='store_true', help='sampler: targets h_l h_l^T, directions = union of the columns of h_l o_l (o_l orthogonal, chosen to overlap) minus one duplicated direction')
ap.add_argument('--plainT', action='store_true', help='positive control: random h_l, T_l = h_l h_l^T, dirs = columns of the h_l (shuffled)')
a = ap.parse_args(); m = a.m
def mm(A, B):
    out = []
    for i in range(m):
        r = 0
        for k in range(m):
            if (A[i] >> k) & 1: r ^= B[k]
        out.append(r)
    return tuple(out)
def tr(A): return tuple(sum(((A[r] >> c) & 1) << r for r in range(m)) for c in range(m))
def rank(A):
    rows = list(A); r = 0
    for c in range(m):
        p = next((i for i in range(r, m) if (rows[i] >> c) & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(m):
            if i != r and (rows[i] >> c) & 1: rows[i] ^= rows[r]
        r += 1
    return r
I = tuple(1 << i for i in range(m))
GLL = [g for g in itertools.product(range(1, 1 << m), repeat=m) if rank(g) == m]
ALLT = sorted({mm(g, tr(g)) for g in GLL})
ORTH = [g for g in GLL if mm(g, tr(g)) == I]
rng = np.random.default_rng(a.seed)
S0 = lms.Slice(m, [I], [1]); NF = S0.F + 1
QVm = np.array([[S0.Qval(Q, x) if x else 0 for x in range(1 << m)] for Q in range(1, NF)])
f = open(a.out, 'a'); t_start = time.time()
for trial in range(a.trials):
    if time.time() - t_start > a.maxtime: break
    if a.share:
        for _ in range(1000):
            hs = [I] + [GLL[rng.integers(len(GLL))] for _ in range(a.W - 1)]
            cols = []
            for h in hs:
                best = None
                for o in ORTH:
                    g = mm(h, o)
                    cs = [sum(((g[i] >> j) & 1) << i for i in range(m)) for j in range(m)]
                    ov = sum(c in cols for c in cs) + rng.random() * 0.5
                    if best is None or ov > best[0]: best = (ov, cs)
                cols += best[1]
            dup = [c for c in set(cols) if cols.count(c) >= 2]
            if dup: break
        T = [mm(h, tr(h)) for h in hs]
        cols.remove(dup[rng.integers(len(dup))])
        cols = cols[:a.M] if len(cols) >= a.M else cols
        a.dirs = ','.join(map(str, cols)); a.M = len(cols)
    elif a.plainT:
        hs = [GLL[rng.integers(len(GLL))] for _ in range(a.W)]
        T = [mm(h, tr(h)) for h in hs]
        a.dirs = ','.join(str(sum(((h[i] >> j) & 1) << i for i in range(m))) for h in hs for j in range(m))
        a.M = a.W * m
    elif a.targets: T = [lms.parse_T(t, m) for t in a.targets.split(',')]
    else: T = [I] + [ALLT[rng.integers(len(ALLT))] for _ in range(a.W - 1)]
    Ot = np.array([[S0.pairQT(Q, t) for t in T] for Q in range(1, NF)]).sum(1)     # odd arrays per form
    found = False
    for _ in range(200):
        if a.dirs: dirs = [int(v) for v in a.dirs.split(',')]
        else:
            dirs = [int(rng.integers(1, 1 << m)) for _ in range(a.M)]
            # hill-climb on the total deletion deficit (replace one move's direction at a time)
            for step in range(400):
                cnt = np.bincount(dirs, minlength=1 << m)
                deficit = np.maximum(Ot - QVm @ cnt, 0).sum()
                if deficit == 0: break
                t = int(rng.integers(a.M)); old = dirs[t]; dirs[t] = int(rng.integers(1, 1 << m))
                cnt2 = np.bincount(dirs, minlength=1 << m)
                if np.maximum(Ot - QVm @ cnt2, 0).sum() > deficit: dirs[t] = old
        cnt = np.bincount(dirs, minlength=1 << m)
        if np.all(QVm @ cnt >= Ot): found = True; break
    if found:
        rng.shuffle(dirs); dirs = [int(v) for v in dirs]
        S = lms.Slice(m, T, dirs); assert S.deletion_ok()
    else:
        rec = dict(W=a.W, M=a.M, T=T, cost=None, note='no deletion-passing direction sequence found')
        f.write(json.dumps(rec) + '\n'); f.flush(); continue
    best = None
    for rs in range(a.restarts):
        x0 = (rng.standard_normal(S.nu) + 1j * rng.standard_normal(S.nu)) * 0.5
        t0 = time.time(); x, cost, it = lms.lm(S, x0, iters=a.iters)
        if best is None or cost < best[0]: best = (cost, x, it, time.time() - t0)
    rec = dict(m=m, W=a.W, M=a.M, T=['.'.join(map(str, t)) for t in T], dirs=dirs, cost=best[0], iters=best[2], secs=round(best[3], 2), seed=a.seed, trial=trial)
    f.write(json.dumps(rec) + '\n'); f.flush()
    if best[0] < a.hit:
        np.save('hit_s_m%dW%dM%d_s%d_t%d.npy' % (m, a.W, a.M, a.seed, trial), best[1]); print('HIT', json.dumps(rec), flush=True)
print('done %d trials in %.0f s' % (trial + 1, time.time() - t_start))
