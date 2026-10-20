# c18_s8sets.py -- labelled certificates in the singular subgame on the hook's triple
# Sigma1 = (2,1^6) + (3,1^5) + (2,2,2,1,1) PLUS one more S8 irreducible that carries the hook's AGL(3,2)-types.
# Generic labels for every dimension pattern (several random subspaces per pattern). Calibration: Sigma1 alone must give
# exactly 2.8 (graph labels), never more (flag theorem).
import numpy as np, itertools, time
from lib514 import *
from libs8 import *
import c17_s8labels as C17
t0 = time.time(); rng = np.random.default_rng(18)
info = C17.info; ph, pK = C17.ph, C17.pK
def set_search(lams, ntrials=8, max_patterns=4000):
    blocks = [info[l][0] for l in lams]; dims = [Y.dim for Y in blocks]; offs = np.cumsum([0] + dims); n = offs[-1]
    Hh = np.zeros((n, n)); KK = np.zeros((n, n))
    for Yr, o in zip(blocks, offs): Hh[o:o + Yr.dim, o:o + Yr.dim] = Yr(ph); KK[o:o + Yr.dim, o:o + Yr.dim] = Yr(pK)
    copies = {}
    for (lam, Yr, o) in zip(lams, blocks, offs):
        for chi, Qs in info[lam][1].items():
            for Q in Qs:
                E = np.zeros((n, Q.shape[1])); E[o:o + Yr.dim] = Q; copies.setdefault(chi, []).append((E, Yr, o, Q))
    aligned = {}
    for chi, cl in copies.items():
        ref = cl[0]; out = [ref[0]]
        for (E, Yr, o, Q) in cl[1:]:
            Qa = C17.align(ref[3], Q, ref[1], Yr); Ea = np.zeros_like(E); Ea[o:o + Yr.dim] = Qa; out.append(Ea)
        aligned[chi] = out
    chis = list(aligned); mults = [len(aligned[c]) for c in chis]
    best = (0, None); count = 0
    pats = list(itertools.product(*[range(m + 1) for m in mults]))
    if len(pats) > max_patterns: pats = [pats[i] for i in rng.choice(len(pats), max_patterns, replace=False)]
    for dims_ in pats:
        gen = any(0 < d < m for d, m in zip(dims_, mults))
        for trial in range(ntrials if gen else 1):
            cols = []
            for c, d, m in zip(chis, dims_, mults):
                if d == 0: continue
                Lsub = np.linalg.qr(rng.standard_normal((m, d)))[0]
                for j in range(d): cols.append(sum(Lsub[i, j] * aligned[c][i] for i in range(m)))
            if not cols: continue
            Yc = np.hstack(cols)
            if Yc.shape[1] in (0, n): continue
            cap = rank(np.hstack([Yc, Hh @ Yc])) - Yc.shape[1]
            if cap == 0: continue
            dem = rank(np.hstack([Yc, KK @ Yc])) - Yc.shape[1]; count += 1
            if dem / cap > best[0] + 1e-9: best = (dem / cap, (dims_, cap, dem))
    return best, mults, count
base = [(2, 1, 1, 1, 1, 1, 1), (3, 1, 1, 1, 1, 1), (2, 2, 2, 1, 1)]
b, m, c = set_search(base)
print('Sigma1 (hook triple):', 'best %.4f' % b[0], b[1], ' type multiplicities', m, ' certificates tried', c, flush=True)
extras = [(3, 3, 2), (5, 3), (4, 3, 1), (3, 2, 2, 1), (4, 2, 1, 1), (5, 2, 1), (3, 2, 1, 1, 1), (7, 1), (6, 1, 1)]
for lam in extras:
    b, m, c = set_search(base + [lam])
    print(f'Sigma1 + {lam}: best {b[0]:.4f} {b[1]}  multiplicities {m}  tried {c}  [{time.time()-t0:.0f}s]', flush=True)
