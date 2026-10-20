"""fano_gate.py (mw2-singer4): ms-fano's Fano-template prefixes (fano7.prefix: PL, FS, FN, FE) fed to the new end test
(tops.py): for each paid element E of the pool, the all-T gate (top augmentations over C[B'], RULES D2) and the Singer
test (all T over C[<eta>], Rado per character).  This is the sparse/fast replacement for ms-fano's joint2.py (which
solved only three fixed 2x2 coefficient matrices, took minutes per prefix and capped on 36/207 systems).
usage: fano_gate.py KIND W k n seed   (block 4; the word has k + W moves)"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../fano-blueprint'))
import numpy as np
import fano7
from s4lib import S4, E_pool, fast_frame_key
from tops import target_space, target_space_full, best_tops, rank_p

class S4on(S4):
    """S4 sharing an existing bm.Game (so prefixes built elsewhere can be tested)"""
    def __init__(self, G, eta_src):
        self.m = G.m; self.n = G.n; self.G = G
        self.rng = random.Random(1); self.FK = {}; self.MONO = {}; self.INV = {}
        self.D1 = G.D1
        self.eta = eta_src.eta
        self.etap = list(eta_src.etap)
        self.etap_id = [G.intern(x) for x in self.etap]
        self.eta_id = self.etap_id[1]

def main():
    kind, W, k, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    t0 = time.time(); rng = random.Random(seed)
    src = S4(4, seed=seed)
    ft = fano7.setup(4, seed)
    S = S4on(ft.G, src)
    agg = {'built': 0, 'fail': 0, 'augrank': {}, 'singer_pass': 0, 'full_cap': 0, 'eta_cap': 0}
    tt = []
    for it in range(n):
        if len(ft.G.MATS) > 40000:
            ft = fano7.setup(4, seed + 1000 + it); S = S4on(ft.G, src)
        gens = fano7.prefix(ft, kind, W, k, rng)
        if gens is None:
            agg['fail'] += 1; continue
        agg['built'] += 1
        t1 = time.time()
        pool = E_pool(S, gens, maxn=12, rng=rng)
        imgs = []
        for E in pool:
            img, st = target_space_full(S, gens, W, E)
            if img is None:
                agg['full_cap'] += 1; continue
            imgs.append(img)
        r = rank_p(np.concatenate(imgs, axis=0)) if imgs else 0
        agg['augrank'][r] = agg['augrank'].get(r, 0) + 1
        # how many prefix generators carry K-adjacent terms (dK = 1)?  record (count, rank)
        nk = sum(1 for g in gens if any(S.G.dist(gid) == 1 for (l, gid) in g))
        key = '%d:%d' % (nk, r)
        agg.setdefault('nK1_rank', {})
        agg['nK1_rank'][key] = agg['nK1_rank'].get(key, 0) + 1
        if r >= W:
            spaces = {}
            for idx, E in enumerate([None] + pool):
                V, st = target_space(S, gens, W, E, S.etap_id)
                if V is None:
                    agg['eta_cap'] += 1; continue
                spaces['E%d' % idx] = V
            combo, diag = best_tops(spaces, W, 7)
            if combo is not None:
                agg['singer_pass'] += 1
                print('SINGER PASS', it, combo, diag, flush=True)
        tt.append(time.time() - t1)
        if it % 20 == 19:
            print('progress', it + 1, agg, 'median s/prefix %.2f' % float(np.median(tt)), flush=True)
    agg['seconds'] = round(time.time() - t0, 1)
    agg['median_s_per_prefix'] = round(float(np.median(tt)), 3) if tt else None
    print(f'fano_gate {kind} W={W} k={k} n={n} seed={seed}: {agg}', flush=True)

if __name__ == '__main__':
    main()
