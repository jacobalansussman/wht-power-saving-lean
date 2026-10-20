#!/usr/bin/env python3
"""d-bit-exits / dirt_join.py: the exact exit obstruction for #589's helpers.

Claim (dirt-column lemma). In a word whose gates are reversible additions at a common frame, with frames that only
grow, a helper h can end at a frame U only if, at every time, every register whose frame is not inside U holds no
component of h's dirt column e_h. (The e_h column restricted to the registers outside U changes only by row operations
among them, since any gate touching one of them lifts all its operands outside U; new registers can join that set but
none can leave; a nonzero vector stays nonzero under invertible row operations and under extension; at the end every
register other than h must hold no e_h.)
So, whatever cleanup follows B (any gates, any order, any frames), h can end below FULL only if
  J_h = join of the end-of-B frames of all registers holding e_h after B (h included) is not FULL.
This script computes J_h for every helper, and also J_h at the time of h's last touch in A/B (frames at that time,
registers holding e_h then): if that is FULL too, no re-timed cleanup of h inside A/B can end it below FULL either
without changing what A and B compute.
usage: python3 -I dirt_join.py WORD.json.gz [out.json]"""
import gzip, json, os, sys
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames
w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R; F = Frames(w)
helpers = range(2 * V, N)
last = {}
seq = [('D', g) for g in w['D']] + [('A', g) for g in w['A']] + [('S', None)] + [('B', g) for g in w['B']]
for i, (sec, g) in enumerate(seq):
    if g is None: continue
    for r in [g[1], *g[2]]: last[r] = i
cur = list(w['start']); val = [1 << i for i in range(N)]
at_last = {}
by_time = {}
for h in helpers: by_time.setdefault(last[h], []).append(h)
for i, (sec, g) in enumerate(seq):
    if g is None:
        snap = [val[q] for q in w['ret']]
        for t, ks in enumerate(w['scatter']):
            for k in ks: val[V + t] ^= snap[k]
    else:
        f, t, ss = g
        for r in [t, *ss]: cur[r] = f
        for s in ss: val[t] ^= val[s]
    for h in by_time.get(i, []):
        holders = [r for r in range(N) if val[r] >> h & 1]
        J = F.join([cur[r] for r in holders])
        at_last[h] = (F.dim[J], len(holders), F.dim[cur[h]])
endJ = {}
for h in helpers:
    holders = [r for r in range(N) if val[r] >> h & 1]
    endJ[h] = (F.dim[F.join([cur[r] for r in holders])], len(holders))
out = dict(end_of_B_join_dim=dict(sorted(Counter(endJ[h][0] for h in helpers).items())),
           end_of_B_holders=dict(sorted(Counter(endJ[h][1] for h in helpers).items())),
           at_last_touch_join_dim=dict(sorted(Counter(at_last[h][0] for h in helpers).items())),
           at_last_touch_holders=dict(sorted(Counter(at_last[h][1] for h in helpers).items())),
           helpers_exit_possible_after_B=sum(1 for h in helpers if endJ[h][0] < 8),
           helpers_with_join_below_FULL_at_last_touch=sum(1 for h in helpers if at_last[h][0] < 8))
print(json.dumps(out, indent=1))
if len(sys.argv) > 2: json.dump(out, open(sys.argv[2], 'w'), indent=1)
