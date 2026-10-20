#!/usr/bin/env python3
"""d-bit-exits / entrances.py: which #589 helpers could START at their first gate frame (an entrance) instead of 0?

A helper that starts at E != 0 cannot be read by D (frame 0), so its target response K_h (the helpers D reads into
targets) must be zero, or be cancelled by donors as in the early-reset kernel (PR435/417). This script lists, for
every helper: its first run (rank of the climb 0 -> first gate frame), whether D reads it, and what an entrance at
the first gate frame would save in rank mass.
Also tests the exit condition exactly (linear algebra over F2 on the post-B z state):
  a helper h can be restored at a hyperplane U >= its end-of-B frame by gates among registers whose end-of-B frames lie
  inside U, with every other z register cleaned at FULL without reading h, only if
  (i) no z register whose end-of-B frame is NOT inside U carries h's dirt column e_h after B, and
  (ii) e_h lies in the F2 span of the post-B values of the registers inside U.
usage: python3 -I entrances.py WORD.json.gz [out.json]"""
import gzip, json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames_lib import Frames

w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R
F = Frames(w)
helpers = range(2 * V, N)
# D reads
Dread = Counter()
for f, t, ss in w['D']:
    for s in ss: Dread[s] += 1
# replay with frames and values (values: z inputs only matter; targets never feed z, H1)
cur = list(w['start']); val = [1 << i for i in range(N)]
first = {}; runs = defaultdict(list)
def gate(f, t, ss):
    for r in [t, *ss]:
        if cur[r] != f:
            assert F.le(cur[r], f)
            runs[r].append((F.dim[cur[r]], F.dim[f], f)); cur[r] = f
    for s in ss: val[t] ^= val[s]
for sec in ('D', 'A'):
    for f, t, ss in w[sec]: gate(f, t, ss)
snap = [val[q] for q in w['ret']]
for t, ks in enumerate(w['scatter']):
    for k in ks: val[V + t] ^= snap[k]
for f, t, ss in w['B']: gate(f, t, ss)
endB = list(cur); postB = list(val)
zregs = [r for r in range(N) if not (V <= r < 2 * V)]

ent = Counter(); ent_zero = Counter(); per = {}
for h in helpers:
    r0 = runs[h][0]                                   # first climb: (0, d1, frame)
    assert r0[0] == 0
    zero_resp = Dread[h] == 0
    ent[r0[1]] += 1
    if zero_resp: ent_zero[r0[1]] += 1
    per[h] = dict(first_dim=r0[1], first_frame=r0[2], D_reads=Dread[h])
out = dict(helpers_read_by_D=sum(1 for h in helpers if Dread[h]), helpers_zero_response=sum(1 for h in helpers if not Dread[h]),
           first_run_rank_all=dict(sorted(ent.items())), first_run_rank_zero_response=dict(sorted(ent_zero.items())),
           entrance_mass_saving_zero_response=sum(d * n for d, n in ent_zero.items()),
           entrance_mass_saving_all_if_responses_cancelled=sum(d * n for d, n in ent.items()))

# exit test (i)+(ii) for X = {h}, U = closure of end-of-B frame of h if dim 7, else every hyperplane above it among
# the end-of-B frames of dim 7 that contain it
def span_contains(vecs, target):
    basis = {}
    for v in vecs:
        x = v
        while x:
            p = x.bit_length() - 1
            if p in basis: x ^= basis[p]
            else: basis[p] = x; break
    x = target
    while x:
        p = x.bit_length() - 1
        if p in basis: x ^= basis[p]
        else: return False
    return True
hyper = sorted({endB[r] for r in zregs if F.dim[endB[r]] == 7})
res = Counter(); okh = []
for h in helpers:
    Us = [endB[h]] if F.dim[endB[h]] == 7 else [U for U in hyper if F.le(endB[h], U)]
    best = None
    for U in Us:
        inside = [r for r in zregs if F.le(endB[r], U)]
        outside = [r for r in zregs if not F.le(endB[r], U)]
        c1 = all(not (postB[r] >> h & 1) for r in outside)
        c2 = span_contains([postB[r] for r in inside], 1 << h)
        key = (c1, c2)
        if best is None or key > best[0]: best = (key, U, len(inside), sum(1 for r in outside if postB[r] >> h & 1))
    if best is None:
        res['no hyperplane above end-of-B frame among dim-7 end frames'] += 1; continue
    (c1, c2), U, ni, nleak = best
    res['(i) %s (ii) %s' % (c1, c2)] += 1
    if c1 and c2: okh.append(h)
    per[h].update(endB_dim=F.dim[endB[h]], exit_i=c1, exit_ii=c2, dirt_leaks_outside_U=nleak)
out['exit_single_helper'] = dict(res)
out['exit_feasible_helpers'] = okh
leak = Counter(per[h].get('dirt_leaks_outside_U', -1) for h in helpers)
out['dirt_leak_count_hist'] = dict(sorted(leak.items()))
print(json.dumps(out, indent=1))
if len(sys.argv) > 2:
    out['per_helper'] = per
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
