#!/usr/bin/env python3
"""an-build: host classes of the slots of a gcert/1 local word under TWO-SIDED slack (s = a 10th local coordinate shared by
all y's; see twosided_local.py for why it is forced).  Own code, stdlib.  usage: python3 -I hostclass.py FILE.gcert.json.gz [label]

Same-run or cross-run hosting on a data array: an X array never contains s (it ends in s^perp); a Y array contains s from
time 0 (its source label).  So a hosted slot either never holds s (X host) or holds s at every gate (Y host).
Pass 1: minimal propagation (s only from the y's and from the totals at the cut).  Pass 2: every slot that never has to
avoid s is put on a Y host (s from its start), the propagation is re-run, and slots are demoted to banks until there is
no violating gate.  Result: per class the number of slots and their residual widths (span = final dim - start dim, +1 for a
slot that must climb s inside the call), which is what the hosting price needs."""
import gzip, json, sys, os
from collections import Counter
path = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(path)
d = json.load(gzip.open(path))
h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
dim = [len(f) for f in d['frames']]
A, B = d['A'], d['B']; GA = A + B; nA = len(A)
tot = [s_ for _, s_, _ in d['ret']]
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
RG = [regs(g) for g in GA]

def propagate(src_from_start):
    """forward s-bit: y's and src_from_start hold s from time 0; totals get s at the cut. returns per-gate bit, final bits"""
    s = [False] * n
    for r in range(v, 2 * v): s[r] = True
    for r in src_from_start: s[r] = True
    at = []
    for i, rs in enumerate(RG):
        if i == nA:
            for t in tot: s[t] = True
        j = any(s[r] for r in rs)
        if j:
            for r in rs: s[r] = True
        at.append(j)
    return at, s

def avoid(x_hosted):
    """backward avoid-bit: x's and X-hosted slots avoid s up to their end"""
    a = [False] * n
    for r in range(v): a[r] = True
    for r in x_hosted: a[r] = True
    at = [False] * len(RG)
    for i in range(len(RG) - 1, -1, -1):
        rs = RG[i]
        j = any(a[r] for r in rs)
        if j:
            for r in rs: a[r] = True
        at[i] = j
    return at

slots = list(range(2 * v, n))
gates_of = {r: [] for r in range(n)}
for i, rs in enumerate(RG):
    for r in rs: gates_of[r].append(i)
G1, _ = propagate([])
A1 = avoid([])
viol1 = sum(1 for i in range(len(RG)) if G1[i] and A1[i])
cls = {}
for r in slots:
    ga = [i for i in gates_of[r] if A1[i]]; gg = [i for i in gates_of[r] if G1[i]]
    holds = bool(gg) or r in tot
    if ga and holds: cls[r] = 'crossing'
    elif ga: cls[r] = 'x'
    elif holds: cls[r] = 'y?'
    else: cls[r] = 'free'
# pass 2: Y-host every 'y?' slot (s from its start); demote to bank until legal
yset = {r for r in slots if cls[r] == 'y?'}
xset = {r for r in slots if cls[r] == 'x'}
it = 0
while True:
    it += 1
    G2, _ = propagate(sorted(yset))
    A2 = avoid(sorted(xset))
    bad = [i for i in range(len(RG)) if G2[i] and A2[i]]
    if not bad: break
    # demote the Y-hosted slots that take part in a violating gate (they will sit on banks and take s when needed)
    culprits = {r for i in bad for r in RG[i] if r in yset}
    if not culprits:
        raise SystemExit('violation without a Y-hosted culprit: %d gates' % len(bad))
    yset -= culprits
    for r in culprits: cls[r] = 'bank-y'
for r in yset: cls[r] = 'y'
# free slots: hostable either side; check they stay s-free under pass 2 (else they become y-type)
G2, s_final = propagate(sorted(yset))
for r in slots:
    if cls[r] == 'free' and any(G2[i] for i in gates_of[r]): cls[r] = 'y-late'
start, final = d['start'], d['final']
def width(r):
    w = dim[final[r]] - dim[start[r]]
    return w + (1 if cls[r] in ('crossing', 'bank-y', 'y-late') else 0)   # these climb s inside the call
C = Counter(cls.values())
Wd = {c: Counter() for c in C}
for r in slots: Wd[cls[r]][width(r)] += 1
res = dict(label=label, file=path, v=v, R=R, N=d.get('N'), cst=d.get('cst'), violations_pass1=viol1, pass2_iterations=it,
           classes=dict(C), widths={c: dict(sorted(w.items())) for c, w in Wd.items()},
           width_sum={c: sum(k * m for k, m in w.items()) for c, w in Wd.items()},
           totals=dict(Counter(cls[t] for t in tot)),
           local_rank=sum(dim[final[r]] - dim[start[r]] for r in range(n)))
print(json.dumps(res))
