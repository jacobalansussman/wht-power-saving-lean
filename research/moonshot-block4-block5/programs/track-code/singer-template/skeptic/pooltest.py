"""skeptic v-mw2-singer4: does the 12-class pool cap (and one_move's randomized legality filter) hide rank?
For FE prefixes at target cost, compute the all-T gate rank with (a) the original pool (maxn=12), (b) the full pool,
(c) the full pool without the one_move filter (every rowsupp-2 E = K g^-1)."""
import sys, random, time
sys.path.insert(0, '.')

import numpy as np
import fano7
from s4lib import S4, E_pool, fast_frame_key
from tops import target_space_full, rank_p
from fano_gate import S4on

def pool_all(S, gens, legal=True):
    G = S.G; out = {}; nrej = 0
    for gen in gens:
        for (l, gid) in gen:
            if G.dist(gid) != 1:
                continue
            E = G.mul(G.K, np.conj(G.MATS[gid]).T)
            k = fast_frame_key(E.T.copy())
            if k in out:
                continue
            if not (G.rowsupp(E) == 2).all():
                continue
            if legal and not G.one_move(E):
                nrej += 1; continue
            out[k] = E
    return list(out.values()), nrej

def rank_of(S, gens, W, pool):
    imgs = []; cap = 0
    for E in pool:
        img, st = target_space_full(S, gens, W, E)
        if img is None:
            cap += 1; continue
        imgs.append(img)
    return (rank_p(np.concatenate(imgs, axis=0)) if imgs else 0), cap

def main():
    W, k, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    src = S4(4, seed=seed); ft = fano7.setup(4, seed); S = S4on(ft.G, src)
    stats = {'built': 0, 'r12': {}, 'rfull': {}, 'rnolegal': {}, 'poolsize': {}, 'onemove_rej': 0, 'diff': 0, 'caps': 0}
    t0 = time.time()
    for it in range(n):
        if len(ft.G.MATS) > 40000:
            ft = fano7.setup(4, seed + 1000 + it); S = S4on(ft.G, src)
        gens = fano7.prefix(ft, 'FE', W, k, rng)
        if gens is None:
            continue
        stats['built'] += 1
        p12 = E_pool(S, gens, maxn=12, rng=random.Random(it))
        pf, nrej = pool_all(S, gens, True)
        pn, _ = pool_all(S, gens, False)
        stats['onemove_rej'] += nrej
        r12, c1 = rank_of(S, gens, W, p12); rf, c2 = rank_of(S, gens, W, pf); rn, c3 = rank_of(S, gens, W, pn)
        stats['caps'] += c1 + c2 + c3
        for key, r in (('r12', r12), ('rfull', rf), ('rnolegal', rn)):
            stats[key][r] = stats[key].get(r, 0) + 1
        stats['poolsize'][len(pf)] = stats['poolsize'].get(len(pf), 0) + 1
        if not (r12 == rf == rn):
            stats['diff'] += 1
            print('DIFF', it, r12, rf, rn, len(p12), len(pf), len(pn), flush=True)
    stats['sec'] = round(time.time() - t0, 1)
    print(f'pooltest FE W={W} k={k} n={n} seed={seed}:', stats, flush=True)

main()
