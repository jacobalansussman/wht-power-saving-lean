#!/usr/bin/env python3
"""v-pair-market-reprice: exact frame-join replay of U8's B + tail with each target lifted to its own plane P_own(t)
right after its rank-2 stop (the partner of that stop is NOT lifted: first read).  Gate frame := F[f] + cur[a] + cur[b].
Options: skip cross-propagation through the 112 M copies (d-mcopy-code mc_groups.json) and/or target-target gates
(those registers then only climb to F[f] + own frame).  Reports targets that leave u_t^perp and the partner class of
the first gate that pushed each one out.
Usage: python3 -I lift_join.py WORD DEMAND BLOCKS MCGROUPS"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
d = json.load(open(sys.argv[2])); U8 = json.load(open(sys.argv[3]))['blocks']
M = {r for g in json.load(open(sys.argv[4])) for r in g['members']}
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']; h = c['h']
blk = {s: i for i, x in enumerate(U8) for s in x}
def B_(gens):
    piv = {}
    for a in gens:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(sorted(piv.values(), reverse=True))
def cont(a, b): return B_(list(a) + list(b)) == B_(b)
perp = {t: B_([x for x in range(1, 1 << h) if bin(x & ports[t]).count('1') % 2 == 0]) for t in range(v)}
dem = {int(t): {int(s): x for s, x in dd.items()} for t, dd in d['demand'].items()}
sub_readers = defaultdict(list)
for t in range(v):
    per = defaultdict(list)
    for s in dem[t]: per[blk[s]].append(s)
    for b, ss in per.items(): sub_readers[tuple(sorted(ss))].append(t)
own = {t: B_([ports[s] for s in ss]) for ss, ts in sub_readers.items() if len(ts) == 4 for t in ts}
print('M copies:', len(M), ' their start dims:', Counter(len(F[c['start'][r]]) for r in M))
def cls(y): return 'target' if v <= y < 2*v else 'M copy' if y in M else 'source' if y < v else 'other helper'
def run(aug, skipM, skipYY):
    cur = [B_(F[f]) for f in c['start']]
    for ph in ['A', 'precondition', 'cut', 'compensation']:
        for e in c.get(ph) or []:
            f = e[-1]
            if f is None: continue
            for r in (e[1:3] if e[0] == 'add' else e[1:2]):
                if cont(cur[r], F[f]): cur[r] = B_(F[f])
    for r, f in zip(c['centers'], c['center_frames']): cur[r] = B_(F[f])
    r2 = set(); first_out = {}
    for ph in ['B', 'tail']:
        for gi, e in enumerate(c.get(ph) or []):
            f = e[-1]
            if f is None: continue
            regs = list(e[1:3] if e[0] == 'add' else e[1:2])
            skip = len(regs) == 2 and ((skipM and (regs[0] in M or regs[1] in M)) or
                                       (skipYY and all(v <= r < 2*v for r in regs)))
            if skip or len(regs) == 1:
                new = {r: B_(list(F[f]) + list(cur[r])) for r in regs}
            else:
                G = B_(list(F[f]) + list(cur[regs[0]]) + list(cur[regs[1]])); new = {r: G for r in regs}
            for r in regs:
                cur[r] = new[r]
                if aug and v <= r < 2*v and len(F[f]) >= 2 and r not in r2:
                    r2.add(r); cur[r] = B_(list(cur[r]) + list(own[r - v]))
                if v <= r < 2*v and r not in first_out and not cont(cur[r], perp[r - v]):
                    other = [y for y in regs if y != r]
                    first_out[r] = (ph, gi, cls(other[0]) if other else 'scale')
    return first_out
for aug, skipM, skipYY in ((False, False, False), (True, False, False), (True, True, False), (True, True, True)):
    fo = run(aug, skipM, skipYY)
    print(f'aug={aug} skipM={skipM} skipYY={skipYY}: targets pushed outside u^perp {len(fo)}; first push by {dict(Counter(x[2] for x in fo.values()))}')
