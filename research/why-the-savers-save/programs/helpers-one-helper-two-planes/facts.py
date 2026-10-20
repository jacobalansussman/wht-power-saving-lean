#!/usr/bin/env python3
"""Small exact facts used in the d-aboveplane report: compensation gates (frame, position, registers touched), block spans
vs target perps, and the arrival-sandwich summary. Usage: python3 -I facts.py FULL-WORD.json.gz survey.json arrival.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); AR = json.load(open(sys.argv[3])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']
ent = {int(q) for q in w.get('start_extra', {})}; low = {int(q) for q in w.get('helper_final', {})}
print('compensation gates', len(comp), 'all at frame []:', all(e[-1] == [] for e in comp), 'program positions', off['compensation'], '..', off['A'] - 1,
      '(all before phase A); readers are targets:', all(v <= e[1] < 2 * v for e in comp), '; read registers entering or low-ended:',
      sum(1 for e in comp if e[2] in ent or e[2] in low))
blocks = {}
for r in S['planes']: blocks.setdefault(r['block'], set()).update(ports[s] for s in r['sources'])
Q = {b: basis(sorted(x)) for b, x in blocks.items()}
print('block spans: ranks', Counter(len(x) for x in Q.values()), '; (target, block) with Q_b inside u_t^perp:', sum(1 for b in Q for u in ports if inperp(Q[b], u)), 'of', len(Q) * v)
print('arrival sandwich dims by class', Counter((r['cls'], r['dimS']) for r in AR['rows']))
