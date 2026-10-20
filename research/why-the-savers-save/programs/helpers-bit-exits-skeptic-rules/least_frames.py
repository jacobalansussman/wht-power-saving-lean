#!/usr/bin/env python3
"""v-bit-exits-rules: (1) are #589's D/A/B gate frames the least frames (join of operand frames)?
(2) recompute the dirt-join obstruction with least frames in D, A, B (a re-framed word can only lower frames).
(3) span test for the helpers whose join is below FULL.
usage: python3 -I least_frames.py WORD.json.gz"""
import gzip, json, os, sys
from collections import Counter
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../helpers-bit-exits')))
from frames_lib import Frames
w = json.load(gzip.open(sys.argv[1], 'rt'))
V, R = w['v'], w['R']; N = 2 * V + R; F = Frames(w)
helpers = range(2 * V, N)
def run(least):
    cur = list(w['start']); val = [1 << i for i in range(N)]; notleast = Counter()
    for sec in ['D', 'A', 'S', 'B']:
        if sec == 'S':
            snap = [val[q] for q in w['ret']]
            for t, ks in enumerate(w['scatter']):
                for k in ks: val[V + t] ^= snap[k]
            continue
        for f, t, ss in w[sec]:
            j = F.join([cur[r] for r in [t, *ss]])
            if j != f: notleast[sec] += 1
            g = j if least else f
            assert F.le(j, f)
            for r in [t, *ss]: cur[r] = g
            for s in ss: val[t] ^= val[s]
    return cur, val, notleast
for least in (False, True):
    cur, val, nl = run(least)
    J = {}
    for h in helpers:
        holders = [r for r in range(N) if val[r] >> h & 1]
        J[h] = F.join([cur[r] for r in holders])
    print('least' if least else 'word', 'gates not at least frame:', dict(nl),
          'J dims:', dict(sorted(Counter(F.dim[J[h]] for h in helpers).items())))
    # span test for sub-FULL joins: e_h in F2-span of values of registers whose frame <= U, for some U >= J of dim 7
    sub = [h for h in helpers if F.dim[J[h]] < 8]
    # candidate U: every existing frame of dim 7 containing J (frames list in word + any closed hyperplane: use word frames)
    hyp = [k for k in range(len(F.mask)) if F.dim[k] == 7]
    def inspan(vecs, target):
        piv = {}
        for v in vecs:
            for b in sorted(piv, reverse=True):
                if v >> b & 1: v ^= piv[b]
            if v: piv[v.bit_length() - 1] = v
        for b in sorted(piv, reverse=True):
            if target >> b & 1: target ^= piv[b]
        return target == 0
    ok = 0
    for h in sub:
        cands = [U for U in hyp if F.le(J[h], U)] or [J[h]]
        good = [U for U in cands if inspan([val[r] for r in range(N) if F.le(cur[r], U)], 1 << h)]
        ok += bool(good)
    print('   sub-FULL joins:', len(sub), 'of which pass span test for some listed U:', ok, '(hyperplane frames listed:', len(hyp), ')')
