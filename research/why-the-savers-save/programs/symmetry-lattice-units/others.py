"""Other Clifford/code systems: (1) complex side, best imaginable (h = clique lower bound, free copies) against E8 P112
at the same s; (2) bit side, the copy charge x a unit with window n needs to tie E8 P111 at the same s."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from ceil_grid import Xf, x_needed
sv = {}
for f in ('survey_e8_bw16_bw16_other_rm14_rm24w4.json', 'survey_bw32_leech.json'):
    sv.update(json.load(open(f)))
S = (3, 5, 7, 9, 11)
out = {}
for name, rec in sv.items():
    if name == 'e8': continue
    w = rec['clique']; n = rec['n']; v = rec['v']
    cx = {s: Xf(s, w, Fr(0)) / Xf(s, 9, Fr(64, 112)) for s in S}
    bit = {s: x_needed(n, Xf(s, 8, Fr(63, 111)), s=s) for s in S}
    out[name] = dict(v=v, n=n, clique=w, complex_best_imaginable_over_E8P112=cx, bit_x_needed_to_tie_E8P111=bit,
                     bit_cst_needed_s9=(bit[9] * v if bit[9] else None))
    print('%-11s v=%6d n=%2d clique>=%3d | complex (h>=clique, x=0)/E8P112 by s: %s | bit x needed by s: %s' % (
        name, v, n, w, {s: round(r, 3) for s, r in cx.items()}, {s: (round(x, 3) if x else 'impossible') for s, x in bit.items()}))
json.dump(out, open('others.json', 'w'), indent=1)
