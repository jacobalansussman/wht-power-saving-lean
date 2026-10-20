"""Kernel numerators (term_le n=3 vs term_le4 n=4) for the tier-2 ceiling tables (what-if: tier 2 needs exit-frame Lean)."""
import sys, os, json, gzip
from collections import Counter
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from route_price import IDLE
from exact_moment import certify
from rate_gap import best_A
spec = json.load(open(sys.argv[1]))
for name, sched in zip(sys.argv[2::2], sys.argv[3::2]):
    s = json.loads(gzip.decompress(open(sched, 'rb').read()))
    c = spec[name]
    prof = Counter({int(k): n for k, n in c['profile'].items()}); prof[1] += 112; prof[8] += 112
    hist = Counter({r: 72 * n for r, n in prof.items()})
    for r, n in IDLE.items(): hist[r] += 112 * n
    hist.update(p for p in s['pads'] if p)
    W = 64 * 112 + len(s['pads'])
    u = dict(m=73, W=W, D=2560, H=dict(hist))
    rate = certify(dict(hist), 73, W); A0 = int(F(rate['root_lower']) * 10 ** 10)
    print(name, 'exact grid', A0, 'n=3 kernel', best_A(u, A0, 3), 'n=4 kernel', best_A(u, A0, 4))
