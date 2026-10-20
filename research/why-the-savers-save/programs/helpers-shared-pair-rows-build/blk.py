#!/usr/bin/env python3
"""pm-build: what the raw U8 word does inside one block (raw numbering). Values: exact, symbolic over all columns,
in the audit's program order. Frames: audit climb rule. Prints, for each register that ever sits at a frame inside the
block span (dims 1..3) or that is read by a block target while carrying block content, its event list."""
import sys, json, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from collections import Counter, defaultdict
from wordlib import load, program, frame_replay, basis, contained, regs, apply, ends
w = load(sys.argv[1]); S8 = [int(x) for x in sys.argv[2].split(',')]
v, R = w['report']['v'], w['report']['R']; n = 2*v+R
P = w['ports']; span8 = basis([P[s] for s in S8]); S8s = set(S8)
prog = program(w); chain, bad = frame_replay(w, prog); print('frame violations', len(bad))
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
def nm(r): return f'{cls(r)}{r}'
# frame of register q at time T
import bisect
def fr(q, T):
    ch = chain[q]; i = 0
    for k, (t, U) in enumerate(ch):
        if t <= T: i = k
        else: break
    return ch[i][1]
def inb(U): return len(U) > 0 and contained(U, span8)
rows = [{k: F(1)} for k in range(n)]
ev = defaultdict(list)
T8 = set(v+s for s in S8)
watch = set()
for q in range(n):
    if any(inb(U) for _, U in chain[q]): watch.add(q)
watch |= T8
def short(row):
    b = {k: x for k, x in row.items() if k in S8s}
    o = Counter(cls(k) for k in row if k not in S8s)
    return '{' + ','.join(f'{k}:{x}' for k, x in sorted(b.items())) + '}' + (f'+{dict(o)}' if o else '')
for T, (ph, e) in enumerate(prog):
    rs = regs(e)
    if any(r in watch for r in rs):
        if e[0] == 'add':
            a, b = rs; q = F(e[3], e[4])
            line = f'{T:6d} {ph:12s} {nm(a)} += {q} * {nm(b)}  fr={len(e[-1]) if e[-1] is not None else None}{"(in)" if e[-1] and inb(basis(e[-1])) else ""}  b={short(rows[b])}'
        else:
            line = f'{T:6d} {ph:12s} {nm(e[1])} *= {F(e[2], e[3])}'
        for r in rs:
            if r in watch: ev[r].append(line)
    apply(e, rows)
start, end = ends(w)
for q in sorted(watch):
    print(f'=== {nm(q)} start dim {len(start[q])} end dim {len(end[q])} chain dims', [len(U) for _, U in chain[q]],
          'inblock', [(t, len(U)) for t, U in chain[q] if inb(U)])
    for l in ev[q]:
        if 'cleanup' in l and len(ev[q]) > 60: continue
        print('   ', l)
