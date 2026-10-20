"""Hub format ceiling X(s, h, k=1, x) (format-limits-corrected skel.py, theorem variant; certified upper ends) for
the candidate units, against the E8 units the hub uses.  Best s from a float scan over odd s in 3..61, then
certified at that s and at the hub's s = 9."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hubcopy'))
from fractions import Fraction as Fr
from skel import X
from ceil_grid import Xf

UNITS = [  # name, side, h, x (cst/v), note
    ('E8 P112 (hub complex)', 'complex', 9, Fr(64, 112), 'hub #523..#1042 unit'),
    ('E8 P111 (hub bit)', 'bit', 8, Fr(63, 111), 'hub #545/#555/#589 unit'),
    ('E8 P120', 'bit', 8, Fr(72, 120), 'full E8'),
    ('BW16 complex, h>=72, copies free (x=0)', 'complex', 72, Fr(0), 'upper bound for every BW16 complex unit with Gram decoder'),
    ('BW16 bit, ovoid E8-block decoder', 'bit', 16, Fr(1296, 2160), 'x = 0.6, constructed'),
    ('BW16 bit, hypothetical 72 full-rank totals', 'bit', 16, Fr(72 * 16, 2160), 'x = 0.533, not constructed'),
    ('BW16 bit, copies free (x=0)', 'bit', 16, Fr(0), 'upper bound for every BW16 bit unit'),
    ('BW32 complex, h>=172, x=0', 'complex', 172, Fr(0), 'upper bound'),
    ('BW32 bit, x=0', 'bit', 32, Fr(0), 'upper bound'),
    ('BW32 bit, x=0.6', 'bit', 32, Fr(3, 5), 'if an E8-block decoder exists'),
    ('Leech complex, h>=357, x=0', 'complex', 357, Fr(0), 'upper bound'),
    ('Leech bit, x=0', 'bit', 24, Fr(0), 'upper bound'),
    ('Leech bit, x=0.6', 'bit', 24, Fr(3, 5), 'if an E8-block decoder exists'),
]
out = []
for name, side, h, x, note in UNITS:
    scan = [(Xf(s, h, x), s) for s in range(3, 62, 2)]
    fb, sb = max(scan)
    rb = X(sb, h, 1, x); r9 = X(9, h, 1, x)
    rec = dict(unit=name, side=side, h=h, x=float(x), note=note, best_s=sb,
               X_best=float(rb['hi']) if rb else None, X_s9=float(r9['hi']) if r9 else None,
               m_s9=int(r9['m']) if r9 else None, D_per_port_s9=float(r9['D']) if r9 else None)
    out.append(rec)
    print('%-48s %-7s h=%3d x=%.4f  X*=%.4e (s=%d)  X(s=9)=%s  m(s=9)=%s' % (name, side, h, float(x), rec['X_best'] or 0, sb,
          '%.4e' % rec['X_s9'] if rec['X_s9'] else 'D<=0', rec['m_s9']), flush=True)
json.dump(out, open('ceilings.json', 'w'), indent=1)
