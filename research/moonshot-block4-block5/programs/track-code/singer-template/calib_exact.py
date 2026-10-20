"""calib_exact.py (mw2-singer4): end-to-end EXACT check of the pipeline on plain-cost hits.
Builds ms-fano FE prefixes (W arrays, k moves) until the all-T gate and the Singer test pass, then solves the two (or W)
tops exactly over Q with s4lib.S4.joint (fixed coefficient matrix alpha from a small menu), rebuilds the whole word from
the solution and checks the FORMAL identity  output_o == K x_o  term by term (exact Fractions, exact group elements).
usage: calib_exact.py W k seed maxhits"""
import sys, os, random, itertools
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../fano-blueprint'))
import fano7
from s4lib import S4, E_pool
from tops import target_space_full, rank_p
from fano_gate import S4on

def alphas(W):
    I = [[({0: 1} if o == j else {}) for j in range(W)] for o in range(W)]
    out = [('identity', I)]
    if W == 2:
        out.append(('halving', [[{0: Fr(1, 2)}, {0: Fr(1, 2)}], [{0: Fr(1, 2)}, {0: Fr(-1, 2)}]]))
        out.append(('one-shared', [[{0: 1}, {0: 1}], [{}, {0: 1}]]))
    return out

def rebuild_and_check(S, gens, W, Es, alpha, sol):
    G = S.G
    tops = [dict() for _ in range(W)]
    contents = [dict() for _ in range(W)]
    outs = [dict() for _ in range(W)]
    for (key, x) in sol:
        if key[0] == 'E':
            _, j, i, b = key
            for t, c in S.lmul(G.MATS[b], gens[i]).items():
                contents[j][t] = contents[j].get(t, 0) + x * c
        else:
            _, o, i, b = key
            for t, c in S.lmul(G.MATS[b], gens[i]).items():
                outs[o][t] = outs[o].get(t, 0) + x * c
    for j in range(W):
        tops[j] = S.lmul(Es[j], {t: c for t, c in contents[j].items() if c != 0})
    for o in range(W):
        acc = dict(outs[o])
        for j in range(W):
            for k, ak in alpha[o][j].items():
                for t, c in S.lmul(S.etap[k], tops[j]).items():
                    acc[t] = acc.get(t, 0) + ak * c
        # raw inputs absorb free terms
        acc = {t: c for t, c in acc.items() if c != 0 and not G.is_monomial(G.MATS[t[1]])}
        if acc != {(o, G.KID): 1}:
            return False, o, len(acc)
    return True, None, sum(len(c) for c in contents)

def main():
    W, k, seed, maxhits = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    src = S4(4, seed=seed); ft = fano7.setup(4, seed); S = S4on(ft.G, src)
    hits = 0; tried = 0
    while hits < maxhits and tried < 2000:
        tried += 1
        gens = fano7.prefix(ft, 'FE', W, k, rng)
        if gens is None:
            continue
        pool = E_pool(S, gens, maxn=12, rng=rng)
        imgs = [target_space_full(S, gens, W, E)[0] for E in pool]
        imgs = [i for i in imgs if i is not None]
        if not imgs or rank_p(np.concatenate(imgs, axis=0)) < W:
            continue
        done = False
        for combo in itertools.product(range(len(pool)), repeat=W):
            Es = [pool[c] for c in combo]
            for name, al in alphas(W):
                r, info = S.joint(gens, W, Es, al)
                if r == 'yes':
                    ok, o, nterms = rebuild_and_check(S, gens, W, Es, al, info['sol'])
                    print(f'prefix {tried}: tops {combo} alpha {name}: exact solve yes; formal identity of the rebuilt '
                          f'{k + W}-move word: {ok} (content terms {nterms})', flush=True)
                    hits += ok; done = True
                    break
            if done:
                break
        if not done:
            print(f'prefix {tried}: gate passed but no exact solution with the alpha menu', flush=True)
    print(f'calib_exact W={W} k={k}: {hits} exact words verified after {tried} prefixes')

if __name__ == '__main__':
    main()
