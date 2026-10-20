"""ballfft.py (w8-b5sparse, resumed): for s1501 (four 8-cycles of coefficients a1-8, a9-16, a17-24, a25-32), rerun the
numerical radius-k LM (ballm's systems) and print the 8-point DFT magnitudes of each cycle of the best alpha (its Fourier
support S = which V4 component the numerical solution lies on).  Hint only.  usage: ballfft.py d4.pkl d3.pkl k starts"""
import sys, pickle, random
import numpy as np
import endsolve as es
import v4gb
from ballm import sub
D4 = pickle.load(open(sys.argv[1], 'rb')); D3 = pickle.load(open(sys.argv[2], 'rb'))
k = int(sys.argv[3]); starts = int(sys.argv[4]); nal = D4['nal']; rng = random.Random(3)
Nss = []
for l in (0, 1):
    cl, sc, cap = D4['sys'][l]; rs, _ = v4gb.leafpeel(cl, sc, None); Nss.append(es.assemble(sub(cl, sc, rs), nal, rng))
    cl, sc, cap = D3['sys'][l]; rs, _ = v4gb.leafpeel(cl, sc, v4gb.ball(cl, sc, k)); Nss.append(es.assemble(sub(cl, sc, rs), nal, rng))
for s in range(starts):
    cost, al = es.lm(Nss, nal, random.Random(100 + s), starts=1, iters=120)
    if al is None:
        continue
    out = []
    for b in range(4):
        v = al[8 * b: 8 * b + 8]
        F = np.abs(np.fft.fft(v)) / max(np.abs(v).max(), 1e-300)
        out.append(' '.join('%.0e' % x for x in F))
    print('start %d residual %.2e min|a| %.2e | DFT/max per cycle: %s' % (s, cost, np.abs(al).min(), ' || '.join(out)), flush=True)
