"""probe_iso.py (g-break, H6) -- complex P112: can a total be charged pairing = star - 1 (drop only to a nonzero
common cut)?  That needs a vector h in span u(R) with u_S . h = 0 for every served S.  Census: for every nonzero h
in F_2^9, the P112 ports with u_S . h = 0 (possible served set), whether h is isotropic, and the largest star-8
read set containing h.  Then the best 'cut decoder' bound: a target S can only be served with discount by totals
whose h is in u_S-perp."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "."))
from common import *
P = list(range(8, 120))
rows = []
for h in range(1, 512):
    Dh = [S for S in P if not dot9(UI[S], h)]
    iso = not dot9(h, h)
    sp = f2dim([UI[S] for S in Dh])
    contains = f2dim([UI[S] for S in Dh] + [h]) == sp     # h in span u(Dh): total reading Dh, serving Dh gets discount
    rows.append((len(Dh), iso, sp, contains, h))
rows.sort(reverse=True)
from collections import Counter
print('h kinds: (|served set in h-perp|, isotropic, span of those labels, h in that span) -> count')
for k, c in sorted(Counter((r[0], r[1], r[2], r[3]) for r in rows).items(), reverse=True):
    print('   ', k, c)
# Is any octad-label h (u_o for o in the removed octad) or u_S for S live special?
oct_h = set(UI[o] for o in range(8))
print('removed-octad labels:', sorted(oct_h))
live_h = set(UI[S] for S in P)
for r in rows[:6]:
    print('  top h %3d: |D| %d iso %s span %d contains %s  (h is a removed-octad label: %s, a live label: %s)' %
          (r[4], r[0], r[1], r[2], r[3], r[4] in oct_h, r[4] in live_h))
# self-served discount family: total reads R = D = P cap h-perp for isotropic h with h in span: pairing = span - 1
disc = [r for r in rows if r[1] and r[3]]
print('isotropic h with h in span of its own served set:', len(disc), ' sizes', Counter((r[0], r[2]) for r in disc))
