"""epscheck.py (v-mw2-absorb4): the eps-family with my own closure (exact E free letters), array 0 (array 1 by symmetry,
checked at the start).  (1) multi-start BFGS over (a2, eps) complex; (2) eps scan (complex eps on circles of radii
1e-4..1e2) with a2 optimised from many starts; (3) coarse grid in a2 at eps = 0.3, 1, 3."""
import sys, pickle, time
import numpy as np
from mycheck import *
from scipy.optimize import minimize
dmp = pickle.load(open(sys.argv[1], 'rb')); m = dmp['m']
isfree = make_free_test('exact', m)
t0 = time.time()
cls = [closure(dmp, l, isfree) for l in range(dmp['W'])]
print("closures", [(len(r), len(c)) for (r, c, e) in cls], f"{time.time()-t0:.0f}s", flush=True)
al = np.array(dmp['alpha'], dtype=complex); fi = dmp['free_idx']; print("free", fi, "alpha", np.round(al, 6))
def tot(v, arrays=(0,)):
    cv = coefs(dmp, v)
    return sum(residual(len(cls[l][0]), len(cls[l][1]), cls[l][2], cv)[0] for l in arrays)
v = al.copy(); v[fi[0]] = 0.7 + 0.2j; v[fi[1]] = 0.3 - 0.1j
print("symmetry check array0 vs array1:", tot(v, (0,)), tot(v, (1,)), flush=True)
rng = np.random.default_rng(3)
def f2(x):
    v = al.copy(); v[fi[0]] = x[0] + 1j * x[1]; v[fi[1]] = x[2] + 1j * x[3]; return tot(v)
print("(1) free starts (array 0 only; total = 2x)", flush=True)
for st in range(int(sys.argv[2])):
    x0 = rng.normal(size=4) * [1, 1, 3, 3]
    r = minimize(f2, x0, method='Nelder-Mead', options={'maxiter': 3000, 'xatol': 1e-9, 'fatol': 1e-12})
    r = minimize(f2, r.x, method='BFGS', options={'maxiter': 200})
    print(f"  start {st}: array-0 res^2 {r.fun:.8f}  a2 {r.x[0]:+.5f}{r.x[1]:+.5f}i  eps {r.x[2]:+.4g}{r.x[3]:+.4g}i", flush=True)
print("(2) eps scan, a2 optimised (10 starts each)", flush=True)
for rad in [1e2, 10, 3, 1, 0.3, 0.1, 1e-2, 1e-3, 1e-4]:
    for ph in [0, 0.5, 1.0, 1.5]:
        eps = rad * np.exp(1j * np.pi * ph)
        def g(x):
            v = al.copy(); v[fi[0]] = x[0] + 1j * x[1]; v[fi[1]] = eps; return tot(v)
        b = min((minimize(g, rng.normal(size=2) * 1.5, method='Nelder-Mead', options={'maxiter': 800}) for _ in range(10)), key=lambda r: r.fun)
        print(f"  eps {rad:g} * exp(i pi {ph}): array-0 res^2 {b.fun:.6f}  a2 {b.x[0]:+.4f}{b.x[1]:+.4f}i", flush=True)
print("(3) a2 grid (re, im in [-3, 3], step 0.25) at eps = 0.3, 1, 3", flush=True)
for eps in [0.3, 1.0, 3.0]:
    best = (9, None)
    for re in np.arange(-3, 3.01, 0.25):
        for im in np.arange(-3, 3.01, 0.25):
            v = al.copy(); v[fi[0]] = re + 1j * im; v[fi[1]] = eps
            c = tot(v)
            if c < best[0]: best = (c, (re, im))
    print(f"  eps {eps}: grid min array-0 res^2 {best[0]:.6f} at a2 {best[1]}", flush=True)
print("done", f"{time.time()-t0:.0f}s")
