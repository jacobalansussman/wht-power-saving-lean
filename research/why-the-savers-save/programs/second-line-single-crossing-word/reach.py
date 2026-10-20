#!/usr/bin/env python3
"""an2-word: time-respecting reachability between the sources x of a gcert/1 local word (own code, stdlib).
Why: with a second free line that DIFFERS from family to family (v_a for family a; X_a never holds v_a), a gate joins equal
frames and frames only grow, so if a time-respecting gate journey leads from x_b's start to x_a, then x_b's start frame
avoids v_a.  If some x_a is reached from EVERY x_b, then no x holds v_a at the start of the call, and the crossing lemma
for the single vector v_a (all y's hold v_a after the scatter, whose connected reader graph forces one common frame on the
readers) needs >= 112 paid crossings of v_a in every GF call: the per-family variant is then no better than shared s
at this gate order.  Gates: phase A, the scatter (copies of the retained totals join the y's), phase B; optionally the
two dirt gates (g1 at time 0: y_T with the slots q, Mt[y_T, q] != 0; g8 at the end: slot q with the registers of its
inverse row).  usage: python3 -I reach.py WORD.json.gz OUT.json"""
import gzip, json, sys
from collections import Counter
d = json.load(gzip.open(sys.argv[1]))
v, R = d['v'], d['R']; n = 2 * v + R
A, B = d['A'], d['B']
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
def reach(with_dirt):
    taint = [0] * n                    # bitmask over the 112 sources: which x starts can have reached this register
    for b in range(v): taint[b] = 1 << b
    seen_x = [0] * v                   # union of the taints x_a carried at its gates
    def gate(rs):
        t = 0
        for r in rs: t |= taint[r]
        for r in rs:
            taint[r] = t
            if r < v: seen_x[r] |= t
    if with_dirt:
        # g1 at time 0: y_T with every slot q of the dirt row (the content of y_T after the main phase depends on q).
        # Taint-wise the y's and slots carry no x yet, so g1 adds no x-to-x journey by itself; it is included for honesty.
        pass
    for g in A: gate(regs(g))
    # scatter: copy k of retained slot s (frame f) is read by the y's of its column; one gate per copy with its readers
    readers = {}
    for T, row in enumerate(d['scat']['table']):
        for k, a, b in row: readers.setdefault(k, []).append(v + T)
    for k, s_, f in d['ret']:
        gate([s_] + readers.get(k, []))
    for g in B: gate(regs(g))
    if with_dirt:
        # g8 at the end: each slot q joins the registers of its inverse row; slots whose final row depends on x join x's.
        # (Exact rows are not needed for reachability: use the forward dependence as a proxy, both are supported on the
        #  same x's for the x-dependent slots of U8X: v-an-build vreplay, end_fwd 68 / end_bwd 354.)  Not modelled here.
        pass
    for a in range(v): seen_x[a] |= taint[a]
    return seen_x
out = {}
for wd in (False,):
    sx = reach(wd)
    sizes = [bin(m).count('1') for m in sx]
    universal = [a for a in range(v) if sizes[a] == v]
    # the best a per-family design could hope for at this order: for each a, the sources that must avoid v_a
    out['reached_from_count_distribution'] = dict(sorted(Counter(sizes).items()))
    out['x_reached_from_all_112'] = len(universal)
    # incoming reach is what matters: x_b must avoid v_a for every a reached from b
    outdeg = [sum(1 for a in range(v) if (sx[a] >> b) & 1) for b in range(v)]
    out['each_x_reaches_count_distribution'] = dict(sorted(Counter(outdeg).items()))
    print('main phase only (A, scatter, B): sources reaching each x (count -> how many x):', out['reached_from_count_distribution'])
    print('  x reached from all 112 sources:', len(universal))
    print('  x_b reaching how many x_a (count -> how many x):', out['each_x_reaches_count_distribution'])
json.dump(out, open(sys.argv[2], 'w'), indent=1)
