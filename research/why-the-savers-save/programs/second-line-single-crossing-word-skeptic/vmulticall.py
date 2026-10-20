#!/usr/bin/env python3
"""v-an2-word: does the crossing lemma's PER-CALL count (>= rank P per call) survive across a SEQUENCE of calls?

Own model and code (stdlib), written independently of an2-word's toybfs.py.
Registers: x (data X array), y (data Y array), and c copies (fresh: start and end at zero).  Status bit per register:
S = its frame holds the shared line s, N = it does not.  Free gate: r_i += lam * r_j only between EQUAL statuses
(a gate joins equal frames).  A status change of a register with nonzero content is a paid crossing (cost 1, the
minimum); a zero register changes status for free.  Contents are vectors over GF(p) in the input basis (x0, y0), so a
solution is the identity on ARBITRARY inputs.
Endpoint contract (two-sided, shared s): x starts N and ends N; y starts S and ends S.
A pair's run is a sequence of lifting calls, e.g. GF, GB, GF (the template's pair 1: stages 0, 1, 2):
   GF: y += x      GB: x -= y.
After call k the data registers must hold exactly the lifting-chain contents (the call is the identity), and copies
must be zero (each call is a closed invocation).
Modes:
  'fixed'  : at every call boundary x is N and y is S (the global frame chain = the old one with s on the Y side) --
             the setting of an2-word's lemma.
  'free'   : statuses at the intermediate call boundaries are free; only the start and the end of the run are fixed.
Minimum total paid crossings by exact 0-1 BFS.
usage: python3 -I vmulticall.py OUT.json"""
import json, sys
from collections import deque


def solve(seq, mode, ncopy, p=3, discard=False):
    n = 2
    def enc(v): return v[0] + p * v[1]
    def dec(c): return [c % p, c // p]
    NV = p * p
    add = [[enc([(dec(i)[0] + dec(j)[0]) % p, (dec(i)[1] + dec(j)[1]) % p]) for j in range(NV)] for i in range(NV)]
    smul = [[enc([(l * dec(i)[0]) % p, (l * dec(i)[1]) % p]) for i in range(NV)] for l in range(p)]
    R = 2 + ncopy
    # lifting chain targets
    x, y = [1, 0], [0, 1]
    targets = []
    for op in seq:
        if op == 'GF': y = [(y[0] + x[0]) % p, (y[1] + x[1]) % p]
        else: x = [(x[0] - y[0]) % p, (x[1] - y[1]) % p]
        targets.append((enc(x), enc(y)))
    start = ((enc([1, 0]), enc([0, 1])) + (0,) * ncopy, 0b10, 0)   # status bit0 = x, bit1 = y ; S = 1
    dist = {start: 0}; dq = deque([start]); best = None
    while dq:
        st = dq.popleft(); c, sb, k = st; d = dist[st]
        if best is not None and d >= best: continue
        if k == len(seq) and (sb & 0b11) == 0b10:
            best = d if best is None else min(best, d); continue
        nxt = []
        # end the current call
        if k < len(seq) and (c[0], c[1]) == targets[k] and (discard or all(v == 0 for v in c[2:])):
            if mode == 'free' or (sb & 0b11) == 0b10:
                # discard: a copy's leftover content is thrown away at the end of the call (it is fresh scratch)
                nxt.append(((c[:2] + (0,) * ncopy, sb, k + 1), 0))
        if k < len(seq):
            for i in range(R):
                si = (sb >> i) & 1
                for j in range(R):
                    if i == j or ((sb >> j) & 1) != si or c[j] == 0: continue
                    for l in range(1, p):
                        nc = list(c); nc[i] = add[c[i]][smul[l][c[j]]]
                        nxt.append(((tuple(nc), sb, k), 0))
        # status changes (allowed any time, also between calls, e.g. an idle climb)
        for i in range(R):
            nxt.append(((c, sb ^ (1 << i), k), 0 if c[i] == 0 else 1))
        for ns, w in nxt:
            if ns not in dist or dist[ns] > d + w:
                dist[ns] = d + w
                (dq.append if w else dq.appendleft)(ns)
    return best, len(dist)


cases = []
seqs = (['GF'], ['GB'], ['GB', 'GF'], ['GF', 'GB'], ['GF', 'GB', 'GF'], ['GF', 'GB', 'GF', 'GB'])
for discard in (False, True):
  for seq in seqs:
    for ncopy in ((1, 2) if not discard else (1,)):
        for mode in ('fixed', 'free'):
            b, ns = solve(seq, mode, ncopy, discard=discard)
            cases.append(dict(seq=seq, copies=ncopy, mode=mode, copies_end='discarded' if discard else 'zero',
                              min_crossings=b, calls=len(seq), states=ns))
            print('%-14s copies %d (%-9s) %-6s  min crossings %s   (calls: %d)   states %d' % ('-'.join(seq), ncopy,
                  'discarded' if discard else 'zero', mode, b, len(seq), ns))
json.dump(cases, open(sys.argv[1], 'w'), indent=1)
