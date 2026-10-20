#!/usr/bin/env python3
"""v-cx-medium / vhost.py (own code, stdlib; python3 -I).  Does hosting make junk cleaning free?  Test of cx-medium's
causal claim ("what sinks the one design that saves is junk cleaning without idle levels"; "hosting is the hub's only
way to make junk cleaning free").  At first order a = D/N; ideal banks are D-neutral, and a hosted helper piece is the
SAME paid piece as on a bank.  So the only difference hosting makes is that the hosted levels REPLACE idle levels of
data arrays (which would otherwise be paid in large idle pieces).  At P = 1 there are no idle levels: so hosting, if it
were possible, would give the same histogram as ideal banks.  Computed here, with cx-medium's own shrink model
(imported, its layout and exact-fill hosting) for P = 2..5, against the same layouts with EVERY helper piece on ideal
banks (no hosting), and the narrow-mass share of the helper pieces in the hosted P = 4 layout (the hub's template)."""
import sys, os, math
from collections import Counter
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../medium-block-networks')))
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
import shrink
from anatomy_ours import price, rank, narrow
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json
cert = json.load(open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X0-checker-histograms.json'))))['U8X0']['hist_by_class']
S_CLS = {int(r): n for r, n in cert['s'].items()}
def nm(H, m): return sum(r * n * math.log(m / r) for r, n in H.items())

for P in (1, 2, 3, 4, 5):
    S, m, ch = shrink.layout(P)
    data = 2 * P * shrink.V
    idle_total = sum(2 * shrink.V * r for _, r, _ in ch)
    T = 200
    # no hosting: every idle child whole (one piece of its rank per instance), all helper levels on ideal banks
    H = Counter({r: S * T * n for r, n in shrink.PROF.items()})
    for _, r, _ in ch: H[r] += 2 * shrink.V * T
    banks = T * S * shrink.LEFT / m
    W = T * data + banks
    D = m * W - rank(dict(H))
    a_nohost = float(price(dict(H), W, m=m)) if D > 0 else 0.0
    hosted = shrink.evaluate(P)
    hs = nm({r: S * n for r, n in S_CLS.items()}, m)          # helper in-call pieces per replica
    print(f"P={P} m={m}: D/replica {D/T:.0f}; hosted (cx-medium exact fill) {hosted['a']:.4e}; all helpers on ideal banks "
          f"{a_nohost:.4e}; helper narrow mass per replica {hs:.0f} = {hs / (hosted['N'] / hosted['T']):.1%} of the hosted N")

# helpers FREE (no helper pieces at all) in the same layouts: the cap that is left when junk cleaning costs nothing
print('helpers free (no helper pieces, idle children whole):')
for P in (1, 2, 4):
    S, m, ch = shrink.layout(P)
    data = 2 * P * shrink.V
    H = Counter({r: S * n for r, n in shrink.PROF.items()})
    for r, n in S_CLS.items(): H[r] -= S * n
    for _, r, _ in ch: H[r] += 2 * shrink.V
    H = {r: n for r, n in H.items() if n}
    W = data; D = m * W - rank(H)
    # hosted levels are now absent: the helper levels were never data levels, so D unchanged
    print(f"  P={P} m={m}: D/replica {D}; price {float(price(H, W, m=m)):.4e}")
