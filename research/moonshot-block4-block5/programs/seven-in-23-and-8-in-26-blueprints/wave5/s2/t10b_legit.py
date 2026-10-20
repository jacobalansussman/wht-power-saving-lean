"""t10b_legit.py (session 2) -- legitimacy of t10's absorbed shared move at the symplectic level:
(1) do the symmetry lifts c_h fix L_Z (are they letters)?  (2) beta_h = m^-1 p_h m c_h^-1 with p_h = RF^-1 pt^h RF:
does sigma(beta_h) fix L_Z (a letter), for every h?  (3) the letters l_X = m^-1 (RF^-1 X RF)^-1 m: do they fix L_Z and
Lam_m?  (4) p is a homomorphism and F1 = m . sum_h psi(h)^-1 beta_h tau_h(y) reproduces t10's terms exactly."""
import sys, os, runpy
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.dirname(here))
from sp8 import *
g = runpy.run_path(os.path.join(here, 't10_build.py'), run_name='t10')
E, RF, mv, pt, grp, basis, y, psi, F1 = g['E'], g['RF'], g['m'], g['pt'], g['grp'], g['basis'], g['y'], g['psi'], g['F1']
P = E.P
def fixes(S, mk): return img_mask(table(S), mk) == mk
LZm = LZ
print('(1) c_h fixes L_Z:', [fixes(E.c[h], LZm) for h in range(7)], '; c\'_h fixes L_Z:', [fixes(E.cp[h], LZm) for h in range(7)])
RFi = inv(RF); mi = inv(mv); m = mv
ph = []; cur = IDM
for h in range(7):
    ph.append(mul(RFi, mul(cur, RF))); cur = mul(pt, cur)
ok_hom = all(mul(ph[a], ph[b]) == ph[(a + b) % 7] for a in range(7) for b in range(7))
beta = [mul(mi, mul(ph[h], mul(m, E.ci[h]))) for h in range(7)]
print('(2) p homomorphism:', ok_hom, '; beta_h fixes L_Z:', [fixes(b, LZm) for b in beta], '; p_h fixes L_Z:', [fixes(p, LZm) for p in ph])
Lam_m = img_mask(table(mi), LZm)
bad = 0
for X in basis[:200]:
    Xs = grp[X]
    lX = mul(mi, mul(inv(mul(RFi, mul(Xs, RF))), m))
    if not (fixes(lX, LZm) and fixes(lX, Lam_m)): bad += 1
print('(3) letters l_X fixing L_Z and Lam_m: bad', bad, 'of', min(200, len(basis)))
# (4) rebuild F1 from the definition
O2 = g['O2']; (l2, R2), = O2.keys()
def rebuild():
    yv = {}
    for c, X in zip(y, basis):
        if abs(c) < 1e-13: continue
        Xs = grp[X]; lX = mul(mi, mul(inv(mul(RFi, mul(Xs, RF))), m))
        u = E.lmul(lX, O2)
        for k, v in u.items(): yv[k] = yv.get(k, 0) + c * v
    acc = {}
    for h in range(7):
        th = E.tau(h, yv)
        bh = E.lmul(beta[h], th)
        for k, v in bh.items(): acc[k] = acc.get(k, 0) + np.conj(E.chars[psi][h]) * v
    return E.lmul(m, {k: v for k, v in acc.items() if abs(v) > 1e-12})
F1b = rebuild()
keys = set(F1) | set(F1b)
diff = max(abs(F1.get(k, 0) - F1b.get(k, 0)) for k in keys)
print('(4) F1 from the definition vs t10 terms: terms', len(F1), len(F1b), ' max coefficient difference %.2e' % diff, flush=True)
