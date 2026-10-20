"""ballm.py (w8-b5sparse, resumed): NUMERICAL hint only (not a verdict): endsolve's LM over alpha on the same systems the
exact one-shot test uses (both outputs' layer-4 closures + their leaf-peeled layer-3 balls of radius k, shared alpha).
A residual near 0 suggests the ball is consistent (so a larger radius is needed); a residual bounded away from 0 over
many starts suggests exclusion.  Controls: s801 k=2 (excluded exactly) must stay away from 0, k=1 must reach 0.
usage: ballm.py sys_d4.pkl sys_d3.pkl k starts iters"""
import sys, pickle, random
import numpy as np
import endsolve as es
import v4gb


def sub(cl, sc, rs):
    rows, cols, ent, kl = cl
    klass, kterms = kl
    R = sorted(rs); rmap = {r: i for i, r in enumerate(R)}
    seedrows = set(rk for c in sc for (rk, cf) in kterms.get(c, []))
    C = sorted(set(c for (r, c, cf) in ent if r in rmap) | set(c for c in klass if any(rk in seedrows for (rk, cf) in kterms.get(c, []))))
    cmap = {c: i for i, c in enumerate(C)}
    ent2 = [(rmap[r], cmap[c], cf) for (r, c, cf) in ent if r in rmap]
    kl2 = ([cmap[c] for c in klass if c in cmap and c in sc], {cmap[c]: [(rk, cf) for (rk, cf) in kterms.get(c, []) if rk in seedrows] for c in sc if c in cmap})
    return ({i: i for i in range(len(R))}, {i: i for i in range(len(C))}, ent2, kl2)


if __name__ == '__main__':
    D4 = pickle.load(open(sys.argv[1], 'rb')); D3 = pickle.load(open(sys.argv[2], 'rb'))
    k = int(sys.argv[3]); starts = int(sys.argv[4]); iters = int(sys.argv[5])
    nal = D4['nal']; rng = random.Random(3)
    Nss = []
    for l in (0, 1):
        cl, sc, cap = D4['sys'][l]
        rs, _ = v4gb.leafpeel(cl, sc, None)
        Nss.append(es.assemble(sub(cl, sc, rs), nal, rng))
        cl, sc, cap = D3['sys'][l]
        rs, _ = v4gb.leafpeel(cl, sc, v4gb.ball(cl, sc, k))
        Nss.append(es.assemble(sub(cl, sc, rs), nal, rng))
    print('systems (rows x cols):', [N[0].shape for N in Nss], flush=True)
    best = []
    for s in range(starts):
        cost, al = es.lm(Nss, nal, random.Random(100 + s), starts=1, iters=iters)
        best.append(cost)
        print('start %d: residual %.3e  min |alpha_k| %.3e' % (s, cost, min(abs(al)) if al is not None else -1), flush=True)
    print('ballm k=%d: best residual %.3e over %d starts' % (k, min(best), starts), flush=True)
