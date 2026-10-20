"""v-mw2-singer4: cross-check the all-T gate (tops.target_space_full) against the independent JOINT end solver
(s4lib.S4.joint: fixed separation matrix alpha, contents solved exactly).  If the gate says rank < W, joint must say no
for EVERY alpha and every choice of tops.  Positive control at plain cost: joint must find yes where the gate passes.
usage: jointcheck.py W k n seed"""
import sys, random, itertools, time
sys.path.insert(0, '.')
import numpy as np
from fractions import Fraction as Fr
import fano7
from s4lib import S4, E_pool
from tops import target_space_full, rank_p
from fano_gate import S4on

def alphas(W):
    out = []
    I = [[({0: 1} if o == j else {}) for j in range(W)] for o in range(W)]
    out.append(('identity', I))
    if W == 2:
        h = Fr(1, 2)
        out.append(('halving', [[{0: h}, {0: h}], [{0: h}, {0: -h}]]))
        out.append(('one-shared', [[{0: 1}, {0: 1}], [{}, {0: 1}]]))
        out.append(('eta-sep', [[{0: 1}, {0: 1}], [{0: 1}, {1: -1}]]))       # x1+x2 / x1 - eta x2 type
        out.append(('eta-plus', [[{0: 1}, {1: 1}], [{1: 1}, {0: 1}]]))
        out.append(('eta-inv', [[{0: h, 1: -h, 2: h}, {}], [{}, {0: 1}]]))     # (1 - eta + eta^2)/2 on one array
    return out

def main():
    W, k, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    src = S4(4, seed=seed); ft = fano7.setup(4, seed); S = S4on(ft.G, src)
    st = {'built': 0, 'gate': {}, 'joint_calls': 0, 'joint': {}, 'contradictions': 0, 'pos_yes': 0}
    t0 = time.time()
    for it in range(n):
        if len(ft.G.MATS) > 40000:
            ft = fano7.setup(4, seed + 1000 + it); S = S4on(ft.G, src)
        gens = fano7.prefix(ft, 'FE', W, k, rng)
        if gens is None:
            continue
        st['built'] += 1
        pool = E_pool(S, gens, maxn=1000)
        imgs = [target_space_full(S, gens, W, E)[0] for E in pool]
        imgs = [i for i in imgs if i is not None]
        r = rank_p(np.concatenate(imgs, axis=0)) if imgs else 0
        st['gate'][r] = st['gate'].get(r, 0) + 1
        if not pool:
            continue
        found = False
        for combo in itertools.combinations_with_replacement(range(len(pool)), W):
            Es = [pool[c] for c in combo]
            for name, al in alphas(W):
                res = S.joint(gens, W, Es, al)
                st['joint_calls'] += 1
                key = res[0] if isinstance(res, tuple) else res
                st['joint'][key] = st['joint'].get(key, 0) + 1
                if key == 'yes':
                    found = True
                    if r < W:
                        st['contradictions'] += 1
                        print('CONTRADICTION', it, combo, name, r, flush=True)
            if found and r >= W:
                break
        if found and r >= W:
            st['pos_yes'] += 1
    st['sec'] = round(time.time() - t0, 1)
    print(f'jointcheck FE W={W} k={k} n={n} seed={seed}:', st, flush=True)

main()
