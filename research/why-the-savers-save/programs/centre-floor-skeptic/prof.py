#!/usr/bin/env python3
"""v-centre / prof.py: own enumeration of the dimension profiles (d_1..d_m), sum <= 63, that pass
  incidence  sum n(d) >= 224 - m,
  NC_k (k <= 6) for the m - k smallest: sum n(d) >= 112 - h(k), h = max E8 lines in rank k (1,3,6,12,20,36),
  NC_7 for the m - 7 smallest: cover Gamma_a (53), N_i (56) by per-dim coverage maxima, or reach 70.
n(d) and the coverage maxima are taken from enum.out (own enumeration)."""
import re, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
txt = open(os.path.join(HERE, 'enum.out')).read()
n = {int(a): int(b) for a, b in re.findall(r'd=(\d)  #W=\d+  n\(d\)=(\d+)', txt)}
cg = {int(a): int(b) for a, b, _ in re.findall(r'cov d=(\d)  Gamma (\d+)\.\.(\d+)', txt)}
cn = {int(a): int(b) for a, b in re.findall(r'cov d=(\d)  Gamma \d+\.\.\d+   N (\d+)', txt)}
h = {0: 0, 1: 1, 2: 3, 3: 6, 4: 12, 5: 20, 6: 36}
print('n', n); print('covG', cg); print('covN', cn)
def ok(ds, m):
    if sum(n[d] for d in ds) < 224 - m: return False
    for k in range(0, 8):
        if m - k <= 0: return False
        J = ds[:m - k]
        if k <= 6:
            if sum(n[d] for d in J) < 112 - h[k]: return False
        else:
            if not (sum(cg[d] for d in J) >= 53 or sum(cn[d] for d in J) >= 56 or sum(n[d] for d in J) >= 70): return False
    return True
res = {}
for m in range(1, 64):
    out = []
    def rec(ds, lo, rem, left):
        if left == 0:
            if ok(ds, m): out.append(tuple(ds))
            return
        for d in range(lo, 10):
            if d * left > rem: break
            rec(ds + [d], d, rem - d, left - 1)
    rec([], 1, 63, m)
    if out: res[m] = out
for m, P in res.items(): print(m, P)
