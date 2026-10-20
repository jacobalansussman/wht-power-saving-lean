#!/usr/bin/env python3
"""pm-build: how U8 delivers each block's content to targets (raw numbering). Symbolic replay over SOURCE columns only
(clean: helpers/targets start 0), audit program order without compensation/cut (as u8_nonplane's demand replay).
For every target read (B phase) of a register carrying source content: the content's block footprint.
Summary: per target, number of content reads, kinds; per block of reader-target-block pairs."""
import sys, json, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from collections import Counter, defaultdict
from wordlib import load, program, regs, apply
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
P = w['ports']
# raw blocks = canonical U8 blocks mapped
m = {int(k): x for k, x in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw2canon_regs.json'))).items()}
inv = {x: k for k, x in m.items()}
CB = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../helpers-shared-pair-rows/u8_blocks.json')))['blocks']
RB = [sorted(inv[s] for s in b) for b in CB]; blk = {s: i for i, b in enumerate(RB) for s in b}
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
rows = [dict() for _ in range(n)]
for s in range(v): rows[s] = {s: F(1)}
prog = [(ph, e) for ph, e in program(w) if ph not in ('compensation', 'cut')]
kinds = Counter(); per_t = defaultdict(list)
for T, (ph, e) in enumerate(prog):
    if e[0] == 'add' and v <= e[1] < 2*v and rows[e[2]]:
        fp = Counter(blk[s] for s in rows[e[2]])
        sig = tuple(sorted(fp.values()))
        kinds[(ph, cls(e[2]), sig)] += 1
        per_t[e[1]].append((T, ph, cls(e[2]), e[2], sig, len(e[-1]) if e[-1] is not None else None))
    apply(e, rows)
print('target content reads by (phase, partner class, block footprint sizes):')
for k, x in sorted(kinds.items(), key=lambda z: -z[1]): print('   ', x, k)
print('reads per target:', Counter(len(x) for x in per_t.values()))
t = int(sys.argv[2]) if len(sys.argv) > 2 else 120
for r in per_t[t]: print('  ', r)
# final target content per block
print('final target content block footprint (first 3 targets):')
for t in range(v, v+3): print('  ', t, sorted(Counter(blk[s] for s in rows[t]).items()))
