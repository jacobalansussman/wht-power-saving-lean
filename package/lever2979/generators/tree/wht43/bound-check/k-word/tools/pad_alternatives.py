"""Check that the single minimal pad beats other pad layouts (same word, same widths) for a price JSON."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from exact_moment import certify
d = json.load(open(sys.argv[1]))
H = Counter({int(k): v for k, v in d['profile'].items()}); W = d['W']; P = int([k for k in d['pad_hist'] if k != '0'][0])
root = lambda h, W: certify(dict(h), 73, W)['root_lower_decimal'][:16]
print('ceil (one pad %d, B=%d):' % (P, d['banks']), root(H, W))
def alt(pads, extra):
    h = H.copy(); h[P] -= 1
    if not h[P]: del h[P]
    for p in pads: h[p] += 1
    return root(h, W + extra)
for pads in [(69, P + 73 - 69), (P + 73 - 57, 57) if P + 73 - 57 <= 69 else (57, P + 16)]:
    print('B+1 banks, pads', pads, alt(pads, 1))
for pads in [(P // 2, P - P // 2), (P - 4, 4)]:
    print('B banks, pads split', pads, alt(pads, 0))
