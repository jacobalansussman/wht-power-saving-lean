"""fs-price: controls of the pricing pipeline (Shea's route price and U8XM's exact price)."""
import sys, os, json, gzip
from collections import Counter
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from route_price import unit, price, widths_from, IDLE
BC = os.environ.get('BC', '.')  # as run: the agents' working folder; see generators/README.txt
words = json.load(open(sys.argv[1]))
out = {}
# 1. Shea: #588 atm3 with 112 source copies, cyclic schedule of dedicated banks (route_price.unit)
c = words['588atm3_shea_route_control']
u = unit(c['profile'], widths_from(c['widths']), copies=True)
p = price(u); out['shea'] = dict(W=u['W'], banks=u['banks'], D=u['D'], **p)
print('Shea control     W', u['W'], 'banks', u['banks'], 'D', u['D'], 'kernel', p['kernel_A'], 'DFT', p['dft'], '(want 12331685)')
# 2. U8XM: k-combine's Lean-checked price file; rebuild its histogram from local profile + idle banks + pads
pf = json.load(open(BC + '/k-combine/out/price/U8XM.json'))
s = json.loads(gzip.decompress(open(BC + '/k-combine/out/price/U8XM.schedule.json.gz', 'rb').read()))
hist = Counter({int(r): 72 * n for r, n in pf['local_profile'].items()})
for r, n in IDLE.items(): hist[r] += n * 112
hist.update(x for x in s['pads'] if x)
W = 64 * 112 + len(s['pads'])
assert {str(k): v for k, v in sorted(hist.items())} == pf['profile'] and W == pf['W']
D = 73 * W - sum(r * n for r, n in hist.items())
p = price(dict(hist=dict(hist), W=W, D=D)); out['u8xm'] = dict(W=W, banks=len(s['pads']), D=D, **p)
print('U8XM control     W', W, 'banks', len(s['pads']), 'D', D, 'kernel', p['kernel_A'], 'DFT', p['dft'], '(want 13868856)')
# 3. C588M (no hosting, m-gen annealed) for the record
pf = json.load(open(BC + '/m-gen/out/price/price-annealed.json'))
hist = Counter({int(r): n for r, n in pf['profile'].items()}) if 'profile' in pf else None
if hist:
    W = pf['W']; D = 73 * W - sum(r * n for r, n in hist.items())
    p = price(dict(hist=dict(hist), W=W, D=D)); out['c588m'] = dict(W=W, D=D, **p)
    print('C588M control    W', W, 'D', D, 'kernel', p['kernel_A'], 'DFT', p['dft'], '(want 13552606)')
json.dump(out, open(sys.argv[2], 'w'), indent=1)
