"""Same-s comparison of the hub ceiling X(s,h,1,x) (float scan of hubcopy/skel.py's profile) between the BW16 units and
the hub's E8 units, for odd s = 3..61.  Prints the worst-case (largest) ratio BW16/E8 over s."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from ceil_grid import Xf, x_needed
S = range(3, 62, 2)
pairs = [('BW16 complex h=72 x=0  vs E8 P112', (72, Fr(0)), (9, Fr(64, 112))),
         ('BW16 bit h=16 x=0.6 (ovoid) vs E8 P111', (16, Fr(3, 5)), (8, Fr(63, 111))),
         ('BW16 bit h=16 x=0.533 (hypothetical) vs E8 P111', (16, Fr(8, 15)), (8, Fr(63, 111)))]
out = {}
for name, (h1, x1), (h0, x0) in pairs:
    rat = [(Xf(s, h1, x1) / Xf(s, h0, x0), s) for s in S]
    worst = max(rat)
    out[name] = dict(max_ratio=worst[0], at_s=worst[1], ratio_s9=dict((s, r) for r, s in rat)[9])
    print('%-50s max over s of X_new/X_E8 = %.3f (s=%d); at s=9: %.3f' % (name, worst[0], worst[1], out[name]['ratio_s9']))
# copy charge a BW16 bit unit would need to tie E8 P111 at each s
need = {s: x_needed(16, Xf(s, 8, Fr(63, 111)), s=s) for s in (3, 5, 7, 9, 11, 13)}
print('x a BW16 bit unit (h=16) needs to tie E8 P111 at the same s:', {s: round(x, 4) if x else None for s, x in need.items()})
out['bw16_bit_x_needed_same_s'] = need
json.dump(out, open('per_s.json', 'w'), indent=1)
