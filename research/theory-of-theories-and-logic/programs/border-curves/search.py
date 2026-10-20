#!/usr/bin/env python3
"""search.py (v-ml-proofcx): my own random-start Levenberg-Marquardt search for a first-order border solution
(A, D, B, C) of a relaxation; prints best residual^2 and parameter sizes.  Usage: search.py FILE NSTARTS TCAP"""
import sys, time, json
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import least_squares
from myfo import parse, fo_system

f, NS, TCAP = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
vs, p, rel = parse(f)
fo = fo_system(vs, rel)
names = sorted({v for P in fo for m in P for v, _ in m}); n = len(names); idx = {v: i for i, v in enumerate(names)}
terms = []  # (eq, coef, [(var, exp)])
for k, P in enumerate(fo):
    for m, c in P.items(): terms.append((k, c, [(idx[v], e) for v, e in m]))
E = len(fo)
def F(z):
    out = np.zeros(E, complex)
    for k, c, m in terms:
        t = c
        for i, e in m: t = t * z[i] ** e
        out[k] += t
    return out
def J(z):
    Jm = np.zeros((E, n), complex)
    for k, c, m in terms:
        for a, (i, e) in enumerate(m):
            t = c * e * z[i] ** (e - 1)
            for b, (i2, e2) in enumerate(m):
                if b != a: t = t * z[i2] ** e2
            Jm[k, i] += t
    return Jm
def fr(x):
    z = x[:n] + 1j * x[n:]; r = F(z); return np.concatenate([r.real, r.imag])
def jr(x):
    z = x[:n] + 1j * x[n:]; Jm = J(z)
    return np.block([[Jm.real, -Jm.imag], [Jm.imag, Jm.real]])
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
t0 = time.time(); best = None
for s in range(NS):
    if time.time() - t0 > TCAP: break
    x0 = rng.normal(size=2 * n)
    try:
        sol = least_squares(fr, x0, jac=jr, method='lm', max_nfev=4000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    except Exception as ex:
        print('start', s, 'error', ex); continue
    r2 = float(np.sum(sol.fun ** 2)); z = sol.x[:n] + 1j * sol.x[n:]
    sz = {pre: float(max(abs(z[i]) for i, v in enumerate(names) if v.startswith(pre))) for pre in ('A_', 'D_', 'B_', 'C_')}
    print(json.dumps({'start': s, 'residual2': r2, 'nfev': sol.nfev, **sz}), flush=True)
    if best is None or r2 < best[0]: best = (r2, z.copy())
    if r2 < 1e-26 and max(sz.values()) < 1e3:
        print('PASS-candidate', s, flush=True)
        json.dump({v: [best[1][i].real, best[1][i].imag] for i, v in enumerate(names)}, open(f.split('/')[-1] + '.sol.json', 'w'))
        break
print('best residual2', best[0] if best else None, 'time %.0f s' % (time.time() - t0))
