# fc-gate route D: the hub moment engine's own enclosures (m-anatomy's exact re-implementation, copy in
# tools/lever/pricer.py of this repository: Mercator logs 224 terms, 24-term exp, 224-bit
# dyadic outward rounding), extended by the second target. Mixed rate fact, divided by W, in the engine's normal form
# (b = 1 - z):   M(b) = sum_r (r n_r)/(m W) (m/r)^b  -  (W - H)/W  -  (H kc)/(m W) (m/kc)^b  <= 0.
# Gives the 1e-24-grid b (largest grid point with M < 0, certified at b and b + 1e-24 as the hub's certify does),
# the largest DFT numerator N with M <= 0 PROVED at b1 = (100 N + 2)/10^12, and N + 1 certified false at (N+1)/10^10.
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), *['..'] * 6, 'tools', 'lever'))
from fractions import Fraction as F
from pricer import log_bounds, exp_bounds, _floor, _ceil, price_b
def M(p, b):
    m, W, H, kc = p['m'], p['W'], p['H'], p['kc']
    lo = hi = F(0)
    for r, n in p['children'].items():
        l, u = log_bounds(F(m, r)); e, f = exp_bounds(b * l, b * u)
        q = F(r * n, m * W); lo += q * e; hi += q * f
    l, u = log_bounds(F(m, kc)); e, f = exp_bounds(b * l, b * u)
    q = F(H * kc, m * W)
    lo -= q * f; hi -= q * e
    c = F(W - H, W); lo -= c; hi -= c
    return _floor(lo), _ceil(hi)
def grid24(p):
    lo, hi = F(0), F(3, 1000)               # M(lo) < 0 < M(hi)
    assert M(p, lo)[1] < 0 and M(p, hi)[0] > 0
    while hi - lo > F(1, 10 ** 24):
        mid = F(((lo + hi) / 2 * 10 ** 26).__floor__(), 10 ** 26)
        a, c = M(p, mid)
        if c < 0: lo = mid
        elif a > 0: hi = mid
        else: raise SystemExit('wide')
    g = F((lo * 10 ** 24).__floor__(), 10 ** 24)
    while M(p, g + F(1, 10 ** 24))[1] < 0: g += F(1, 10 ** 24)
    assert M(p, g)[1] < 0 < M(p, g + F(1, 10 ** 24))[0], 'not certified'
    return g
def scanN(p, g):
    N = int(g * 10 ** 10)
    while not (M(p, F(100 * N + 2, 10 ** 12))[1] <= 0): N -= 1
    while M(p, F(100 * (N + 1) + 2, 10 ** 12))[1] <= 0: N += 1
    return N
base = {1: 8041680, 2: 2306880, 3: 729648, 4: 353160, 5: 69120, 8: 25920}
W = 322560
def mk(extra, H, kc):
    ch = dict(base)
    for r, n in extra.items(): ch[r] = ch.get(r, 0) + n
    return dict(children=ch, W=W, H=H, kc=kc, m=73)
P = {'full_B': mk({24: 114624}, 80640, 25), 'half_B': mk({24: 80640, 48: 57312}, 80640, 49),
     'full_pub': mk({24: 161280}, 57312, 25)}
out = {}
u = price_b({**base, 48: 137952}, 73, W)
out['uniform_hub_engine_b'] = '%.24f' % float(0) if False else str(u.numerator * 10 ** 24 // u.denominator)
for k, p in P.items():
    g = grid24(p); N = scanN(p, g)
    out[k] = dict(b_grid24='0.' + str(g.numerator * 10 ** 24 // g.denominator).rjust(24, '0'), N=N,
                  N_proved=M(p, F(100 * N + 2, 10 ** 12))[1] <= 0,
                  N1_false_at_z1=M(p, F(100 * (N + 1) + 2, 10 ** 12))[0] > 0,
                  N1_false_at_grid=M(p, F(N + 1, 10 ** 10))[0] > 0)
    print(k, out[k], flush=True)
print('uniform (hub engine price_b, x1e-24):', out['uniform_hub_engine_b'])
json.dump(out, open(sys.argv[1], 'w'), indent=1)
