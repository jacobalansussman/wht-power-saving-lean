# vc_1434.py -- exact psi (mod 2 primes) and phi (sum of |entries|^2, dyadic, exact in float64) for orbit-union cuts
import sys, time
import numpy as np
from fractions import Fraction as F
from vclib import *
t0 = time.time()
orb, lab, gens = orbits_gl42()
for b in [int(x) for x in sys.argv[1].split(',')]:
    inS = np.array([(b >> int(lab[k])) & 1 == 1 for k in range(NT)])
    out = {}
    for j in (1, 3):
        Bf = blocks(j, 'f')
        phi = 0.0
        for rows, G in Bf.values():
            m = inS[rows]
            if m.all() or (~m).all(): continue
            phi += float((np.abs(G[np.ix_(~m, m)]) ** 2).sum())
        ps = []
        for p in (998244353, 1000000009):
            ps.append(psi_modp(blocks(j, 'p', p), inS, p))
        out[j] = (ps, phi)
        print(f"bits {b} orbits {[k for k in range(12) if b >> k & 1]}: j={j} psi mod p = {ps}, phi = {phi!r} = {F(phi).limit_denominator(1 << 30)} ({time.time()-t0:.0f}s)", flush=True)
    p1, f1 = out[1][0][0], F(out[1][1]).limit_denominator(1 << 30); p3, f3 = out[3][0][0], F(out[3][1]).limit_denominator(1 << 30)
    print(f"  => ratio3 = {F(p3, p1)} = {p3/p1:.6f}, spec3 = phi3/phi1 = {f3/f1}, s1 = {f1/p1} = {float(f1/p1):.6f}, s3 = {f3/p3} = {float(f3/p3):.6f}", flush=True)
