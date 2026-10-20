#!/usr/bin/env python3
"""v-an2-word: the bank case over a whole run, own model and code (stdlib), exact 0-1 BFS over GF(2).
Registers: x (X array, N -> N), y (Y array, S -> S), a bank b (N at the start, S at the end: it climbs s once in its
life, its content restored), and one copy (fresh scratch, discarded at the end of each call).  Free gates only between
equal statuses; a status change of a nonzero register costs 1; zero registers change status freely.  A run of lifting
calls (the template's pair-1 run GF, GB, GF; over GF(2) GB is x += y).  After each call the data registers hold the
lifting-chain contents and the bank holds b.
Modes: 'fixed' (x N, y S at every call boundary; the bank may climb at any time) and 'free' (data statuses free at the
intermediate boundaries).  The whole-run information bound predicts: fixed = calls + 1, free = 2 + 1 (the bank's own
climb never carries an x channel).
usage: python3 -I vbank.py OUT.json"""
import json, sys
from collections import deque


def solve(seq, mode):
    # contents: 3-bit vectors over (x0, y0, b0); registers 0 x, 1 y, 2 bank, 3 copy; status bit 1 = S
    X0, Y0, B0 = 1, 2, 4
    x, y = X0, Y0; targets = []
    for op in seq:
        if op == 'GF': y ^= x
        else: x ^= y
        targets.append((x, y))
    start = ((X0, Y0, B0, 0), 0b0010, 0)
    dist = {start: 0}; dq = deque([start]); best = None
    while dq:
        st = dq.popleft(); c, sb, k = st; d = dist[st]
        if best is not None and d >= best: continue
        if k == len(seq) and (sb & 0b0111) == 0b0110 and c[2] == B0:
            best = d if best is None else min(best, d); continue
        nxt = []
        if k < len(seq) and (c[0], c[1]) == targets[k] and c[2] == B0:
            if mode == 'free' or (sb & 0b11) == 0b10:
                nxt.append(((c[:3] + (0,), sb, k + 1), 0))
        if k < len(seq):
            for i in range(4):
                for j in range(4):
                    if i != j and ((sb >> i) & 1) == ((sb >> j) & 1) and c[j]:
                        nc = list(c); nc[i] ^= c[j]; nxt.append(((tuple(nc), sb, k), 0))
        for i in range(4):
            nxt.append(((c, sb ^ (1 << i), k), 1 if c[i] else 0))
        for ns, w in nxt:
            if ns not in dist or dist[ns] > d + w:
                dist[ns] = d + w; (dq.append if w else dq.appendleft)(ns)
    return best, len(dist)


res = []
for seq in (['GF'], ['GB'], ['GF', 'GB'], ['GF', 'GB', 'GF']):
    for mode in ('fixed', 'free'):
        b, ns = solve(seq, mode)
        res.append(dict(seq=seq, mode=mode, min_crossings_including_bank=b, states=ns))
        print('%-10s %-6s min crossings (bank included) %s   states %d' % ('-'.join(seq), mode, b, ns))
json.dump(res, open(sys.argv[1], 'w'), indent=1)
