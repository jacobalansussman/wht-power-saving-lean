#!/usr/bin/env python3
"""g-recombine: value map of the FRAGMENTATION lever (same helpers, fewer stops) against the COUNT lever, on #942.
WHAT-IF prices, not constructions.  Per-helper run lists come from the exact U8 replay (stops.py's climb rule);
x, y and centre runs are kept; global = 333 x local + residual.  When helper width is unchanged the residual is T37's
{48:11086, 1:6}; when width falls, the freed width is re-packed smallest-capacity-first (optimistic) as in
u-helpers-hub/price_removals.py.
Usage: python3 -I price_frag.py CANONICAL-WORD.json.gz"""
import sys, os, json, gzip
sys.path.insert(0, os.environ.get('PRICER_DIR', '.'))
from pricer import price_b, kappa_from, C_PR589, dec
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
B = [basis(f) for f in F]; dims = [len(b) for b in B]
def inside(a, b):
    bb = B[b]
    for x in B[a]:
        for z in bb: x = min(x, x ^ z)
        if x: return False
    return True
cur = list(c['start']); runs = defaultdict(list)
def climb(r, f, ph):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(cur[r], f): return
    runs[r].append((dims[cur[r]], dims[f], ph)); cur[r] = f
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f, 'centres')
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1], ph)
        if e[0] == 'add': climb(e[2], e[-1], ph)
for r in range(n): climb(r, c['final'][r], 'final')
cen = set(c['centers'])
def local(hrun_fn):
    prof = Counter(); hw = 0
    for r in range(n):
        rr = runs[r]
        if r >= 2*v:
            rr2 = hrun_fn(r, rr); hw += sum(b - a for a, b, _ in rr2)
        else: rr2 = rr
        for a, b, _ in rr2: prof[b - a] += 1
    return prof, hw
P942 = {1: 848490, 2: 112554, 3: 148518, 4: 53946, 5: 23976, 8: 2664, 48: 11086}; M, W = 73, 33152
b0 = price_b(P942, M, W)
def glob(prof, hw):
    G = Counter({r: 333*x for r, x in prof.items()})
    G[8] += 333*8                                   # the 8 rank-8 centre copy children (hub profile)
    hosted = 37*9*hw
    for cap, nb in ((8, 672*37), (16, 448*37), (32, 448*37), (48, 448*37)):
        full = min(nb, hosted//cap); hosted -= full*cap; left = nb - full
        if left and hosted: G[cap-hosted] += 1; left -= 1; hosted = 0
        G[cap] += left
    assert hosted == 0
    return {r: x for r, x in G.items() if x}
ident = lambda r, rr: rr
p0, hw0 = local(ident)
G0 = glob(p0, hw0)
print('control local', dict(sorted(p0.items())), 'helper width', hw0)
print('control global (optimistic re-pack of T37 bins):', dict(sorted(G0.items())))
def merge_if(pred):
    def f(r, rr):
        out = []
        for a, b, ph in rr:
            if out and out[-1][1] == a and pred(out[-1], (a, b, ph)): out[-1] = (out[-1][0], b, ph)
            else: out.append((a, b, ph))
        return out
    return f
cases = {
 'one run per helper': merge_if(lambda p, q: True),
 'tails >= dim 6 in one run (B/tail/cleanup)': merge_if(lambda p, q: p[0] >= 6),
 'B-phase + cleanup tails from dim 3 up merged': merge_if(lambda p, q: p[0] >= 3 and q[2] in ('B', 'tail', 'cleanup')),
 'factory 0->3 in one run (phase A)': merge_if(lambda p, q: q[1] <= 3 and q[2] == 'A' and p[2] == 'A'),
 'cleanup steps only merged': merge_if(lambda p, q: q[2] == 'cleanup' and p[2] == 'cleanup'),
 'cold helpers enter at first gate frame (drop first run if 0->1)': (lambda r, rr: rr[1:] if rr and rr[0][0] == 0 and rr[0][1] == 1 else rr),
 'pairs of rank-1 runs fused (1+1 -> 2)': None,
}
def pairs(r, rr):
    out = []; i = 0
    while i < len(rr):
        a, b, ph = rr[i]
        if b - a == 1 and i+1 < len(rr) and rr[i+1][0] == b and rr[i+1][1] - rr[i+1][0] == 1:
            out.append((a, rr[i+1][1], ph)); i += 2
        else: out.append((a, b, ph)); i += 1
    return out
cases['pairs of rank-1 runs fused (1+1 -> 2)'] = pairs
_seen = [0]
def pairfac(r, rr):
    # pair factory (mutable sources pre-combine x_a += x_b on their line): a helper whose chain starts 0->1->2 in
    # phase A picks the pair up in one rank-2 run; applied to every other such helper (6 of 12 planes per block)
    if len(rr) >= 2 and rr[0][:2] == (0, 1) and rr[1][:2] == (1, 2) and rr[0][2] == rr[1][2] == 'A':
        _seen[0] += 1
        if _seen[0] % 2: return [(0, 2, 'A')] + rr[2:]
    return rr
cases['pair factory: 0->1->2 fused for half the factory chains'] = pairfac
cases['hyperplane ends: drop the cleanup 8->9 climb (end at the dim-8 join)'] = (lambda r, rr: rr[:-1] if rr and rr[-1][2] == 'cleanup' and rr[-1][:2] == (8, 9) else rr)
_c32 = [0]
def c32(r, rr):
    # stand-in: 32 of the cold helpers lose their first 0->1 run (the promo_scan singletons have dim-1 first frames)
    if rr and rr[0][:2] == (0, 1) and _c32[0] < 32: _c32[0] += 1; return rr[1:]
    return rr
cases['32 cold helpers enter at dim 1 (stand-in for the promo_scan singletons)'] = c32
res = {}
_cold1 = set()
for r in range(2*v, n):
    pass
for name, fn in cases.items():
    p, hw = local(fn); G = glob(p, hw)
    b = price_b(G, M, W); k = kappa_from(b, C_PR589)[0]
    _seen[0] = 0; _c32[0] = 0
    hr = sum(len(fn(r, runs[r])) for r in range(2*v, n)); _seen[0] = 0; _c32[0] = 0
    res[name] = dict(helper_runs=hr,
                     helper_width=hw, b=dec(b, 13), gain_pct=float((b/b0 - 1)*100), kappa=dec(k, 16))
    print(name, json.dumps(res[name]))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'price_frag.json'), 'w'), indent=1)
