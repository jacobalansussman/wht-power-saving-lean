#!/usr/bin/env python3
"""d-bit-exits / price_whatif.py: what a width-7 bit helper would be worth on #589 (a WHAT-IF: no admissible construction).

Scenario: N of #589's 439 helpers lose one rank-1 run and one unit of width in every call (an exit at a hyperplane, or
an entrance at a line; the local rank mass drops by N per call, the deficit 10272 is unchanged, so the idle remainder
grows by 288 N). Three hosting models for the freed width:
  r7    freed units paid as rank-7 children (pessimistic bracket of the brief);
  r42   freed units paid as rank-42 children (optimistic bracket of the brief);
  win   an explicit feasible re-pack in PR #442's real nine-stage windows (windows9.py): width-7 residuals of stages
        2, 4, 6 move into the rank-7 windows 1, 3, 5 (which host nothing in #589; a width-7 residual fills one exactly,
        so its rank-7 child disappears), taken from the rank-42/28 windows first so whole copies empty; everything else
        stays in the frozen #589 selection with 7 in place of 8. This is a feasible hosting, not the MILP optimum.
Prices with pricer.py (exact moment engine, surcharge on, 1e-24 grid) and kappa for b = #942's b and 1.10 x that b.
usage: python3 -I price_whatif.py"""
import json, os, sys
from collections import Counter
from fractions import Fraction as Fr
sys.path.insert(0, __import__('os').environ.get('PRICER_DIR', '.'))  # folder holding pricer.py (it ships with the lever package's tools)
from pricer import price_b, kappa_from, dec
HERE = os.path.dirname(os.path.abspath(__file__))
win = {int(k): v for k, v in json.load(open(os.path.join(HERE, 'windows9.json'))).items()}
sel = json.load(open(os.path.join(HERE, 'data/selection9-p111.json')))['bit']
LOC = Counter({1: 1323, 2: 1342, 3: 79, 4: 44, 5: 102, 6: 11, 7: 19}); T = 288; REPL = 32; m, W = 64, 28416
REM = Counter({2: 13920, 4: 14208, 6: 14208, 7: 21312, 42: 288})
B942 = Fr(1664727258740678855772, 10**24)

def glob(loc, rem):
    G = Counter()
    for r, n in loc.items(): G[r] += T * n
    for r, n in rem.items(): G[r] += n
    G = {r: n for r, n in G.items() if n}
    assert m * W - sum(r * n for r, n in G.items()) == 10272
    return G

def rem_win(N):
    """explicit re-pack, per replica; returns the global remainder (x32)"""
    # residuals per (window, stage) from the frozen selection; widths: N of each stage are width 7
    hosted = Counter()
    for i, w_, s, n in sel['stages']: hosted[i, s] += n
    sevens_left = {s: N for s in range(9)}
    content = {i: [] for i in win}                       # widths placed in window i
    moved = Counter()
    for s, target in ((2, 1), (4, 3), (6, 5)):
        k = min(N, 222)
        moved[s] = k; sevens_left[s] -= k
        content[target] += [7] * k
        # take them out of the largest-rank windows hosting stage s
        for i in sorted([i for (i, ss) in hosted if ss == s], key=lambda i: -win[i][1]):
            t = min(k, hosted[i, s]); hosted[i, s] -= t; k -= t
        assert k == 0
    for (i, s), n in sorted(hosted.items()):
        a = min(n, sevens_left[s]); sevens_left[s] -= a
        content[i] += [7] * a + [8] * (n - a)
    assert all(v == 0 for v in sevens_left.values()), sevens_left
    rem = Counter()
    for i, (p, r, st) in win.items():
        items = sorted(content[i], reverse=True); copies = 222
        per = max(len(x[1]) for x in sel['patterns'] if x[0] == i) if any(x[0] == i for x in sel['patterns']) else 0
        per = per or 1
        # fill copies with up to 'per' residuals (as the frozen pattern) while the sum fits; empty copies stay whole
        used = 0
        while items:
            row = []
            while items and len(row) < max(per, 1) and sum(row) + items[0] <= r: row.append(items.pop(0))
            if not row:
                raise AssertionError(('does not fit', i, items[:3]))
            used += 1
            if r - sum(row): rem[r - sum(row)] += REPL
        assert used <= copies, (i, used)
        rem[r] += REPL * (copies - used)
    rem[7] += 0
    for i, r in [(9, 0)]: pass
    return Counter({r: n for r, n in rem.items() if n})

assert rem_win(0) == REM, rem_win(0)
c0 = price_b(glob(LOC, REM), m, W, True)
assert c0 == Fr(1733427039589844234749, 10**24)
rows = []
for N in (1, 12, 50, 100, 220, 439):
    loc = Counter(LOC); loc[1] -= N
    freed = T * N
    r7 = Counter(REM); r7[7] += freed // 7; r7[freed % 7] += 1 if freed % 7 else 0
    r42 = Counter(REM); r42[42] += freed // 42; r42[freed % 42] += 1 if freed % 42 else 0
    rw = rem_win(N)
    row = dict(N=N)
    for name, rem in (('r7', r7), ('r42', r42), ('win', rw)):
        c = price_b(glob(loc, rem), m, W, True)
        k1 = kappa_from(B942, c)[0]; k2 = kappa_from(B942 * Fr(11, 10), c)[0]
        row[name] = dict(c=dec(c, 13), gain_pct=float((c / c0 - 1) * 100), kappa_b942=dec(k1, 16), kappa_b942x1_10=dec(k2, 16))
    rows.append(row)
    print(json.dumps(row))
k0 = kappa_from(B942, c0)[0]; k0b = kappa_from(B942 * Fr(11, 10), c0)[0]
out = dict(control=dict(c=dec(c0, 25), kappa_b942=dec(k0, 16), kappa_b942x1_10=dec(k0b, 16)), rows=rows,
           note='WHAT-IF prices: no admissible construction exists on #589 (see README).')
print(json.dumps(out['control']))
json.dump(out, open(os.path.join(HERE, 'price_whatif.json'), 'w'), indent=1)
