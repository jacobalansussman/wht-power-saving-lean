"""x-pack side check: the bank-packing ceiling for d-lean-route's tier-2 what-if (#799 / #869-870-890 exit frames, E1,
with source copies).  Widths from d-lean-route's annealed schedule; contents by pack_ilp (ceiling banks, one pad);
table validated with Shea's checks; priced with route_price.price (same rate generator).  What-if: tier 2 needs
new Lean (exit-frame helpers) before any of this is a Lean figure.
Usage: tier2_ceiling.py <d-lean-route sa3 schedule> <spec name in d-lean-route work/variants2.json> <variants2.json> <out>"""
import sys, os, json, gzip
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pack_ilp
from pack_exact import validate
from route_price import price, IDLE
sched, name, spec, out = sys.argv[1:5]
s = json.loads(gzip.decompress(open(sched, 'rb').read()))
res = s['widths']
rows, pads, contents = pack_ilp.build(res)
validate(rows, pads, res)
c = json.load(open(spec))[name]
prof = Counter({int(k): n for k, n in c['profile'].items()}); prof[1] += 112; prof[8] += 112
hist = Counter({r: 72 * n for r, n in prof.items()})
for r, n in IDLE.items(): hist[r] += 112 * n
hist.update(p for p in pads if p)
W = 64 * 112 + len(pads); D = 73 * W - sum(r * n for r, n in hist.items())
p = price(dict(hist=dict(hist), W=W, D=D))
print(name, 'banks', len(pads), 'pads', sorted(Counter(x for x in pads if x).items()), 'W', W, 'D', D, p)
open(out, 'wb').write(gzip.compress(json.dumps(dict(rows=rows, pads=pads, widths=res), separators=(',', ':')).encode(), mtime=0))
