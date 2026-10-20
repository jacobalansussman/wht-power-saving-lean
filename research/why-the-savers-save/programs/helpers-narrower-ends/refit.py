#!/usr/bin/env python3
"""d-hpe-multigate: exact frame refitting by coordinate descent (frame annotations only).
For every A, B and early-cleanup (tail) gate, with all other frames fixed, the gate frame F may be any subspace with
join(previous frames of its operands) <= F <= meet(next frames of its operands); only dim F changes the paid runs.
Pick the dimension that minimises the first-order move cost sum r*ln(73/r) of the affected runs; repeat to a fixed point.
Centre interface frames, starts and ends are fixed.  Used after build_ends.build() (which fixes the new ends).
"""
import math, json
from collections import Counter
from cones_lib import basis, contained, registers
from build_ends import meet
M = 73
def wcost(d): return 0.0 if d == 0 else d*math.log(M/d)
def refit(w2, rounds=6):
    v, h, R = w2['report']['v'], w2['report']['h'], w2['report']['R']; n = 2*v+R
    full = basis([1 << i for i in range(h)])
    end = [full]*v+[basis([x for x in range(1, 1 << h) if bin(x & p).count('1') % 2 == 0]) for p in w2['ports']]+[full]*R
    for q, U in w2['helper_final'].items(): end[int(q)] = basis(U)
    start = [(p,) for p in w2['ports']]+[()]*(v+R)
    for q, U in w2['start_extra'].items(): start[int(q)] = basis(U)
    ev = []   # [kind, idx, regs, frame]
    for i, e in enumerate(w2['A']): ev.append(['A', i, list(registers(e)), basis(e[-1])])
    for k, q in enumerate(w2['centers']): ev.append(['centre', k, [q], basis(w2['center_frames'][k])])
    for i, e in enumerate(w2['B']): ev.append(['B', i, list(registers(e)), basis(e[-1])])
    for k, (i, U) in enumerate(w2['early_cleanup']): ev.append(['tail', k, list(registers(w2['cleanup'][i])), basis(U)])
    seq = [[] for _ in range(n)]   # event ids per register
    for j, x in enumerate(ev):
        for r in x[2]: seq[r].append(j)
    pos = {}
    for r in range(n):
        for t, j in enumerate(seq[r]): pos[(r, j)] = t
    def prevf(r, j):
        t = pos[(r, j)]; return start[r] if t == 0 else ev[seq[r][t-1]][3]
    def nextf(r, j):
        t = pos[(r, j)]; return end[r] if t == len(seq[r])-1 else ev[seq[r][t+1]][3]
    moved = 0
    for _ in range(rounds):
        step = 0
        for j, x in enumerate(ev):
            if x[0] == 'centre': continue
            rs = x[2]; P = [prevf(r, j) for r in rs]; N = [nextf(r, j) for r in rs]
            L = basis([y for p in P for y in p]); U = N[0]
            for z in N[1:]: U = meet(U, z)
            assert contained(L, x[3]) and contained(x[3], U)
            def cost(k): return sum(wcost(k-len(p))+wcost(len(nn)-k) for p, nn in zip(P, N))
            k0 = len(x[3]); best = min(range(len(L), len(U)+1), key=lambda k: (cost(k) - 1e-12*(k == k0), abs(k-k0)))
            if cost(best) < cost(k0)-1e-12:
                pool = list(meet(x[3], U))+list(U); F = L
                for y in pool:
                    if len(F) == best: break
                    G = basis(list(F)+[y])
                    if len(G) > len(F): F = G
                assert len(F) == best and contained(L, F) and contained(F, U)
                x[3] = F; step += 1
        moved += step
        if not step: break
    for x in ev:
        if x[0] in ('A', 'B'): w2[x[0]][x[1]][-1] = list(x[3])
        elif x[0] == 'tail': w2['early_cleanup'][x[1]][1] = list(x[3])
    return moved
def recount(w2):
    v, h, R = w2['report']['v'], w2['report']['h'], w2['report']['R']; n = 2*v+R
    full = basis([1 << i for i in range(h)])
    end = [full]*v+[basis([x for x in range(1, 1 << h) if bin(x & p).count('1') % 2 == 0]) for p in w2['ports']]+[full]*R
    for q, U in w2['helper_final'].items(): end[int(q)] = basis(U)
    start = [(p,) for p in w2['ports']]+[()]*(v+R)
    for q, U in w2['start_extra'].items(): start[int(q)] = basis(U)
    cur = start.copy(); paid = Counter()
    def climb(q, U):
        U = basis(U); assert contained(cur[q], U); d = len(U)-len(cur[q])
        if d: paid[d] += 1
        cur[q] = U
    for section in ('A', 'B'):
        if section == 'B':
            lev = {int(t): basis(U) for t, U in w2.get('target_levels', {}).items()}
            for q, U in zip(w2['centers'], w2['center_frames']):
                U = basis(U); assert cur[q] == U
                served = {lev.get(e[1]-v, ()) for e in w2['cut'] if e[2] == q}
                for T in sorted(served, key=lambda x: -len(x)):
                    if T != U: paid[len(U)-len(T)] += 1
                    U = T
        for e in w2[section]:
            for q in registers(e): climb(q, e[-1])
    for i, U in w2['early_cleanup']:
        for q in registers(w2['cleanup'][i]): climb(q, U)
    for q in range(n): climb(q, end[q])
    inv = Counter(len(end[q])-len(start[q]) for q in range(2*v, n))
    w2['report']['local_profile'] = {str(k): x for k, x in sorted(paid.items())}
    w2['report']['width_inventory'] = {str(k): x for k, x in sorted(inv.items())}
    return dict(sorted(paid.items())), sum(len(end[q])-len(start[q]) for q in range(2*v, n))
def S(prof): return sum(x*wcost(r) for r, x in prof.items())
