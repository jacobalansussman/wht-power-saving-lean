#!/usr/bin/env python3
"""g-small: exact scalar check of 'pre-merge half reads'. Ambient F2^5. Packet plane P = <e1,e2,e3>, labels
u1=e1, u2=e2, u3=e3, u4=e1+e2+e3 (any proper subspace of P holds <= 2 labels). Packet p = x1+2x2+3x3+5x4.
Targets t1, t2 with caps G1 = P+<e4>, G2 = P+<e5> (incomparable; each needs p; final frames u_t^perp = G_i here).
Helpers h1, h2 start at frame 0 with ARBITRARY dirt g1, g2 (symbols); targets carry dirt d1, d2.
Claim tested: P is visited by exactly two registers (h1, h2), yet both targets receive p exactly, helpers end
restored, sources unchanged, with every gate at a common frame and all chains monotone. t2 never sits at P: it
reads h1's half at the line <e1,e2> BEFORE the merge and h2's half at G2 after it (jump-over L -> G2).
Usage: python3 -I halfread_toy.py"""
from fractions import Fraction as Fr
def span(vs):
    bb = []
    for y in vs:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return frozenset(sorted(bb))
def contains(G, H):  # H subset of G
    for x in H:
        y = x
        for z in sorted(G, reverse=True): y = min(y, y ^ z)
        if y: return False
    return True
def closure(b):
    s = {0}
    for x in b: s |= {y ^ x for y in s}
    return frozenset(s)
e = [1 << i for i in range(5)]
lab = {'x1': e[0], 'x2': e[1], 'x3': e[2], 'x4': e[0] ^ e[1] ^ e[2]}
regs = ['x1', 'x2', 'x3', 'x4', 'y1', 'y2', 'h1', 'h2']
val = {r: {('x' if r[0] == 'x' else 'd' if r[0] == 'y' else 'g', r): Fr(1)} for r in regs}
frame = {r: span([lab[r]]) if r in lab else span([]) for r in regs}
visits = {r: [closure(frame[r])] for r in regs}
def add(a, b, k, G):
    G = span(G)
    assert contains(G, frame[a]) and contains(G, frame[b]), (a, b)
    for r in (a, b):
        if closure(G) != closure(frame[r]): frame[r] = G; visits[r].append(closure(G))
    for key, x in val[b].items():
        y = val[a].get(key, 0) + k*x
        if y: val[a][key] = y
        else: val[a].pop(key, None)
P = [e[0], e[1], e[2]]; L12 = [e[0], e[1]]; L34 = [e[2], lab['x4']]
G1 = P + [e[3]]; G2 = P + [e[4]]
c = {'x1': 1, 'x2': 2, 'x3': 3, 'x4': 5}
# differential pre-reads of helper dirt (targets and helpers at frame 0: no climb)
add('y1', 'h1', -1, []); add('y1', 'h2', -1, []); add('y2', 'h1', -1, []); add('y2', 'h2', -1, [])
# halves
add('h1', 'x1', 1, [e[0]]); add('h1', 'x2', 2, L12)           # h1 = g1 + x1 + 2x2 at <e1,e2>
add('h2', 'x3', 3, [e[2]]); add('h2', 'x4', 5, L34)           # h2 = g2 + 3x3 + 5x4 at <e3,u4>
add('y2', 'h1', 1, L12)                                        # t2 reads the first half at the line (pre-merge)
add('h1', 'h2', 1, P)                                          # the merge: h1 = g1 + g2 + p at P (visitors h1, h2)
add('y1', 'h1', 1, G1)                                         # t1 reads the full copy at G1
add('y2', 'h2', 1, G2)                                         # t2 reads the second half at G2 (jumps over P)
# restore helpers (cleanup at large frames; sources climb, as in the ladder)
FULL = [e[i] for i in range(5)]
add('h1', 'h2', -1, FULL); add('h1', 'x1', -1, FULL); add('h1', 'x2', -2, FULL)
add('h2', 'x3', -3, FULL); add('h2', 'x4', -5, FULL)
want_p = {('x', s): Fr(k) for s, k in c.items()}
ok_t1 = val['y1'] == {**{('d', 'y1'): Fr(1)}, **want_p}
ok_t2 = val['y2'] == {**{('d', 'y2'): Fr(1)}, **want_p}
ok_h = all(val[h] == {('g', h): Fr(1)} for h in ('h1', 'h2')) and all(val[x] == {('x', x): Fr(1)} for x in lab)
Pc = closure(span(P))
vis_P = [r for r in regs if Pc in visits[r]]
print('t1 gets d1 + p:', ok_t1, ' t2 gets d2 + p:', ok_t2, ' helpers and sources restored:', ok_h)
print('registers that ever sit at P:', vis_P)
print('t2 frame chain dims:', [len(span([x for x in f if x][:])) and len(f).bit_length() - 1 for f in visits['y2']])
print('t1 frame chain dims:', [len(f).bit_length() - 1 for f in visits['y1']], ' t2:', [len(f).bit_length() - 1 for f in visits['y2']])
print('target final frames inside their caps:', contains(span(G1), frame['y1']) and contains(span(G2), frame['y2']))
