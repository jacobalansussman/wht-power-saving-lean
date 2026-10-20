"""t_full.py: the all-T relaxation gate (target_space_full) on the plain control, then on Singer HW prefixes."""
import sys, random, time
from h7 import *
from tops import target_space_full, rank_p
W, M, n, seed, tmpl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
S = S4(4, seed=seed); rng = random.Random(seed); B = Builder(S, rng)
t0 = time.time(); hist = {}; caps = 0
for it in range(n):
    if len(S.G.MATS) > 40000:
        S = S4(4, seed=seed + 1000 + it); B = Builder(S, rng)
    gens = build(B, tmpl, W, M)
    if gens is None or isinstance(gens, str):
        hist[str(gens)] = hist.get(str(gens), 0) + 1
        continue
    pool = E_pool(S, gens, maxn=12, rng=rng)
    imgs = []
    for E in pool:
        img, st = target_space_full(S, gens, W, E)
        if img is None:
            caps += 1; continue
        imgs.append(img)
    tot = np.concatenate(imgs, axis=0) if imgs else np.zeros((0, W))
    r = rank_p(tot)
    hist[r] = hist.get(r, 0) + 1
    if it < 3 or r >= W:
        print(it, 'pool', len(pool), 'aug rank over all tops', r, [int(rank_p(i)) for i in imgs], flush=True)
print(f'{tmpl} W={W} M={M} n={n}: top-augmentation rank histogram {hist}, caps {caps}, {time.time()-t0:.0f}s')
