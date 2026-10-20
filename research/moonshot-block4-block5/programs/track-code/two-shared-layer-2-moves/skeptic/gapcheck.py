"""v-w5-twoshared gapcheck: replays w5-twoshared's cross / sigma cores (their core builders, same seeds) and tests two
hypotheses their frame-level scripts rely on but do not check:
 (G1) cross cores (rd2X / rd2Xb / esearchX): no two layer-3 terms of v on one array share a frame (if they did, letters
      of eps_e could cancel each other's images and the per-term production condition would be too strict);
 (G2) the production source u_k -> t_e with the SAME merge k (k2 == k1, skipped in prod_cond) is only possible for letters
      e' with S_e'(L_D) = L_D; the skip is safe iff no J-condition letter e' (the set cj) lies in Stab(L_D).
Also prints the layer of every main term q_kj.
Usage: gapcheck.py cross|sig1|sig0 ncores seed ..."""
import sys, os, random
W5 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'their')
sys.path.insert(0, W5)
import numpy as np
import sp
from esearch import maps, basis, img_L0, GL, SYM
from sweep import core
from sweep_sig import core_sig
D = sp.D
LD = img_L0(D)
STLD = maps(basis(LD), LD)          # e' in Stab(L0) with S_e'(L_D) = L_D
print('|Stab(L0) cap Stab(L_D)| =', int(STLD.sum()), ' (expect 20,643,840 / 30 = 688,128)', flush=True)

def single(S):
    out = np.zeros((len(GL), 1024), dtype=bool)
    if sp.lag(S) != sp.L0: return out
    # brute: find index by matching (slow path avoided: use symp_from over GL rows with matching P)
    P = tuple(S[i] & 15 for i in range(4))
    i = GL.index(P) if P in GL else None
    for b in range(1024):
        if sp.symp_from(list(GL[i]), SYM[b]) == tuple(S):
            out[i, b] = True; break
    return out

variant = sys.argv[1]; ncores = int(sys.argv[2]); seeds = [int(x) for x in sys.argv[3:]]
TOT = {'G1': 0, 'G2': 0, 'cores': 0}
for seed in seeds:
    rng = random.Random(seed)
    for c in range(ncores):
        if variant == 'cross':
            cr = core(rng, 'cross')
            if cr is None: print('core failed'); continue
            R0, R1, a11, a20, tg = cr
        else:
            cr = core_sig(rng, variant)
            if cr is None: print('core failed'); continue
            R0, R1, a11, a20, tg, o = cr
        TOT['cores'] += 1
        b = {0: R0['b'][0], 1: R1['b'][0]}; s = {0: R0['s'], 1: R1['s']}
        g = {1: R0['g'][0], 2: R1['g'][0]}
        a = {(1, 0): R0['a'][0], (2, 1): R1['a'][0]}
        if a20 is not None: a[(2, 0)] = a20[0]
        if a11 is not None: a[(1, 1)] = a11[0]
        p = {kj: sp.mv(a[kj], s[kj[1]]) for kj in a}
        q = {kj: sp.mv(g[kj[0]], p[kj]) for kj in a}
        mu = {0: R0['mu'][0], 1: R1['mu'][0]}; route = {0: (1, 0), 1: (2, 1)}
        lay = {kj: sp.lam(q[kj]) for kj in q}
        lay3 = [kj for kj in q if lay[kj] == 3]
        shared = []
        for j in (0, 1):
            fr = [sp.lag(q[kj]) for kj in lay3 if kj[1] == j]
            if len(fr) != len(set(fr)): shared.append(j)
        g1 = bool(shared) and variant == 'cross'
        TOT['G1'] += g1
        line = f"{variant} seed {seed} core {c}: tags {tg}, lambda(q) {dict(sorted(lay.items()))}, shared layer-3 frames on arrays {shared}"
        g2hits = []
        for e in (0, 1):
            y = sp.comp(D, mu[e])
            for kj in q:
                if kj == route[e]: continue
                J = sp.comp(y, q[kj]); L = sp.lam(J)
                if L != 2: continue
                arr = kj[1]
                pfr = {sp.lag(p[(k2, arr)]) for k2 in (1, 2) if (k2, arr) in p}
                if sp.lag(J) in pfr: continue
                cj = maps(basis(img_L0(D)), img_L0(sp.comp(J, sp.inv(s[arr]))))
                for l in (1, 2):
                    cj |= maps(basis(img_L0(sp.comp(D, sp.comp(g[l], D)))), img_L0(J))
                for kj2 in q:
                    if kj2[1] == arr and lay[kj2] == 2:
                        cj |= single(sp.comp(J, sp.inv(q[kj2])))
                n = int((cj & STLD).sum())
                if n: g2hits.append((e, kj, n, int(cj.sum())))
        TOT['G2'] += bool(g2hits)
        print(line + f"; G1 {'HIT' if g1 else 'ok'}; G2 (top, junk, |cj cap Stab(L_D)|, |cj|): {g2hits if g2hits else 'none'}", flush=True)
print('TOTAL', TOT, flush=True)
