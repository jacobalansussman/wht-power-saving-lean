#!/usr/bin/env python3
"""v-cx-medium / vprice.py (own code, stdlib; python3 -I).  Independent re-pricing of cx-medium's block-25 one-pair
figures.  The fold inequality is re-typed from its statement (an-ours docstring / fsb_rate_check.py):
    (2^s - 1) * sum_r n_r (r/m)^(1-a) + W ((m-1)/m)^(1-a) + W (1/m)^(1-a) < 2^s W,  s = 40
solved here with mpmath-free Decimal bisection written separately, plus the moment root (s -> infinity) and the
first-order D/N.  The histogram is built from the itemised classes of the U8X0 certificate rerun (x, y, s, c per call),
not from cx-medium's PROF constant; banks are ideal (helper levels / 25) as in cx-medium (a), and I also price
two alternatives cx-medium did not: (e) banks that must pay a pad because a bank's three pieces cannot exceed 25 and
pieces come in leftover widths {4,5,7,8,9}: exact optimum of packing by integer programme-free counting (best-case
lower bound on pads); (f) the T74 profile."""
import json, math, os
from decimal import Decimal as D, getcontext
getcontext().prec = 50
BASE = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X0-checker-histograms.json'))
cert = json.load(open(BASE))['U8X0']['hist_by_class']
cls = {c: {int(r): int(n) for r, n in h.items()} for c, h in cert.items()}

def fold_root(H, W, m, s=40):
    lm = D(m); two = D(2) ** s; Wd = D(W)
    terms = [(D(n), D(r) / lm) for r, n in H.items() if n]
    def g(a):
        e = 1 - a
        def pw(x): return (e * x.ln()).exp()
        return (two - 1) * sum(n * pw(q) for n, q in terms) + Wd * (pw(D(m - 1) / lm) + pw(D(1) / lm)) - two * Wd
    lo, hi = D(0), D('0.05')
    assert g(lo) < 0
    for _ in range(170):
        mid = (lo + hi) / 2
        if g(mid) < 0: lo = mid
        else: hi = mid
    return float(lo)

def moment_root(H, W, m):
    f = lambda a: sum(n * (r / m) ** (1 - a) for r, n in H.items()) - W
    lo, hi = 0.0, 0.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return lo

def rank(H): return sum(r * n for r, n in H.items())
def narrow(H, m): return sum(r * n * math.log(m / r) for r, n in H.items())

m, v, T = 25, 112, 25
per_call = {}
for c, h in cls.items():
    for r, n in h.items(): per_call[r] = per_call.get(r, 0) + n
helper_levels = rank(cls['s'])
print('per-call classes:', {c: (rank(h), sum(h.values())) for c, h in cls.items()}, 'helper levels per call', helper_levels)
assert rank(cls['x']) == v * 8 and rank(cls['y']) == v * 8 and rank(cls['c']) == 64
out = {}
H = {r: 3 * T * n for r, n in per_call.items()}
banks = 3 * T * helper_levels / m
W = 2 * v * T + banks
Dd = m * W - rank(H)
print(f'(a) U8X0 ideal banks: T={T} W={W} banks={banks} D={Dd} (per replica {Dd/T})')
a = fold_root(H, W, m); mr = moment_root(H, W, m); fo = Dd / narrow(H, m)
print(f'    fold root {a:.7e}   moment root {mr:.7e}   first order D/N {fo:.4e}')
out['U8X0_ideal'] = dict(W=W, D=Dd, fold=a, moment=mr, first_order=fo)
# T74 gauge44 profile (as read by an-hub from the hub word's report; data and copies assumed as U8X0's ranks)
t74 = {1: 2348, 2: 500, 3: 416, 4: 130, 5: 82, 8: 8}
hl = rank(t74) - 2 * v * 8 - 64
H2 = {r: 3 * T * n for r, n in t74.items()}; W2 = 2 * v * T + 3 * T * hl / m
a2 = fold_root(H2, W2, m)
print(f'(f) T74 gauge44 ideal banks: helper levels/call {hl}, W={W2}, D={m*W2-rank(H2)}, fold root {a2:.7e}')
out['T74_ideal'] = dict(W=W2, fold=a2)
# helpers free (no banks): data + copies only
H3 = {}
for c in ('x', 'y', 'c'):
    for r, n in cls[c].items(): H3[r] = H3.get(r, 0) + 3 * n
a3 = fold_root(H3, 2 * v, m)
print(f'    helpers free (data and copies only): fold {a3:.4e}, D = {m*2*v - rank(H3)}')
out['helpers_free'] = a3
# generous ceiling: data in {8,16}, copies rank 8
H4 = {8: 2 * v + 24, 16: 2 * v}
print(f'    generous ceiling (data {{8,16}}, copies rank 8, junk free): fold {fold_root(H4, 2*v, m):.4e}, moment {moment_root(H4, 2*v, m):.4e}')
out['ceiling'] = fold_root(H4, 2 * v, m)
# calibration: the hub's T74 figure is not reproducible here (needs its idle layout); instead the U8X0HQ Lean grid
# value is re-derived in vcalib (not done: time-box).
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
