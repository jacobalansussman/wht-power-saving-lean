"""foofo.py -- e3-bpexact: the group-algebra construction in the order F O O F O (7 in 23: F1 O1 O2 F2 O3).
F1 = uniform psi-twisted shared START: m_s . sum_h psi(h)^-1 c_h^-1 tau_h(x_0)  (absorbed in every direction; at the frame
     level it lives in character psi only, Theorem R).
O1 = m1 . (x_0 + b F1)  (two routes),  O2 = m2 . O1,
F2 = group-algebra shared move on O2 (same letters/twist/design as the survivor's F1), character psi too,
O3 = the top (free).  Character psi (F2's own) is to be served through the F1 route: F1 (d1) -> O1 (d2) -> O2 (d3, next
to the kernel) -> O3; the other six through F2's group-algebra pattern.  m_s and b are searched so that the F1 route's
O2 term lies at input distance 3 next to the kernel cell.  Frame-sum relaxation of the lifted skeleton, ALL 7 characters.
argv: seed psi nBm"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
from relax import Eng
t0 = time.time()
seed, psi, nBm = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 1); nrng = np.random.default_rng(seed + 2)
ys = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'ys_seed{seed}_psi{psi}.npy'))
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < nBm:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
P = L.P
# search m_s (a move D(z) letter) and the letter b: cell(m2 m1 b m_s) at distance 3, adjacent to K
found = None
srng = random.Random(seed + 11)
for it in range(40000):
    ms = mm(G.D[srng.randrange(1, 16)], srng.choice(L.letters)[0])
    b = srng.choice(L.letters)[0]
    w = mm(L.m2, L.m1, b, ms)
    c = L.cell(w)
    if P.dZ[c] == 3 and P.D[c, L.iK] == 1 and P.dZ[L.cell(mm(L.m1, b, ms))] == 2:
        found = (ms, b); break
assert found is not None, 'no F1 route found'
ms, b = found
print(f'[foofo] F1 route found after {it + 1} tries: cells d1 {L.cell(ms)}, d2 {L.cell(mm(L.m1, b, ms))}, d3 {L.cell(mm(L.m2, L.m1, b, ms))}', flush=True)
ps = np.conj(L.chars[psi])
F1 = {}
for h in range(7): L.fs_add(F1, h, mm(ms, L.cp[h]), ps[h])          # m_s c_h^-1 tau_h(x_0) = [m_s c'_h] x_h
# absorbed: tau_g(F1) = psi(g) c_g F1
worst = 0
for g in range(7):
    lhs = L.tau(g, F1); rhs = L.lmul(L.c[g], F1); rhs = {k: [L.chars[psi][g] * v[0], v[1]] for k, v in rhs.items()}
    worst = max(worst, L.fs_diff(lhs, rhs))
X0 = {}; L.fs_add(X0, 0, I16, 1.0)
O1 = L.lmul(L.m1, L.clean({**X0, **L.lmul(b, F1)}))
O2 = L.lmul(L.m2, O1)
# F2: group-algebra shared move on O2 (y = sum_b y_b b O2), psi-twisted: m . sum_h psi(h)^-1 beta_h tau_h(y)
y = {}
for yy, bb in zip(ys, Bm):
    for k, (c, w) in O2.items(): L.fs_add(y, k[0], mm(bb, w), yy * c)
acc = {}
for h in range(7):
    th = L.lmul(L.beta[h], L.tau(h, y))
    for k, (c, w) in th.items(): L.fs_add(acc, k[0], w, ps[h] * c)
F2 = L.clean(L.lmul(L.m, L.clean(acc)))
print(f'[foofo] F1 {len(F1)} terms (absorption {worst:.0e}), O1 {len(O1)}, O2 {len(O2)}, F2 {len(F2)} terms '
      f'(absorption {L.check_absorbed(F2):.0e})  {time.time() - t0:.0f}s', flush=True)
E = Eng('Z7', P=L.P, seed=seed)
def conv(u):
    out = {}
    for (l, kk), (c, w) in u.items():
        k2 = (l, inv(sigma(w))); out[k2] = out.get(k2, 0) + c
    return {k: v for k, v in out.items() if abs(v) > 1e-12}
O1c, O2c, F1c, F2c = conv(O1), conv(O2), conv(F1), conv(F2)
for label, shared in (('F O O F O (both)', [F1c, F2c]), ('without F1 route (F2 only)', [F2c])):
    ok, info = E.decide([E.inp(0), O1c, O2c], shared, E=[])
    print(f'[foofo] FRAME-SUM, ALL 7 characters, {label}: passes = {ok}; mu_ok {info.get("mu_ok")}', flush=True)
# control: O1 without the F1 route
O1n = L.lmul(L.m1, X0); O2n = L.lmul(L.m2, O1n)
ok, info = E.decide([E.inp(0), conv(O1n), conv(O2n)], [F1c, F2c], E=[])
print(f'[foofo] control, O1 not reading F1: passes = {ok}; mu_ok {info.get("mu_ok")}', flush=True)
print('done %.0fs' % (time.time() - t0))
