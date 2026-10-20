#!/usr/bin/env python3
"""an2-word: exhaustive check of the crossing lemma on tiny words (own code, stdlib).

Model (the frame game reduced to the one bit that matters): every register has a status, S (its frame holds the extra
line s) or N (it does not).  A free gate may only combine registers of EQUAL status:  r_i += lam * r_j.  A status change
of a register is a paid frame change (a climb N->S or a drop S->N, at least one move), except on a register that holds
zero (free: a fresh register may start at any frame).  x registers are always N, y registers always S (two-sided
endpoint frames); a slot (hosted helper) must start and end at given statuses with its content restored; a copy is a
fresh register (starts at 0, ends with anything).  Exact search (0-1 BFS) over all contents in GF(p)^n.
Goal GF:  (x, y, q) -> (x, y + P x, q);  goal GB: (x, y, q) -> (x - P^T y, y, q).
The lemma predicts: minimum number of paid status changes = rank P + (number of slots whose status differs from start
to end), and the slots' own changes (their 'natural' climbs) never pay for the x -> y channels.
usage: python3 -I toybfs.py OUT.json"""
import itertools, json, sys
from collections import deque


def run(case):
    p, nx, ny = case['p'], case['nx'], case['ny']
    slots = case['slots']                  # list of (start_status, end_status)
    ncopy = case['copies']
    P = case['P']; mode = case['mode']
    nq = len(slots)
    n = nx + ny + nq                       # inputs
    R = n + ncopy
    N = p ** n
    def enc(v):
        c = 0
        for t in reversed(v): c = c * p + t
        return c
    def dec(c):
        v = []
        for _ in range(n): v.append(c % p); c //= p
        return v
    vecs = [dec(c) for c in range(N)]
    add = [[enc([(a + b) % p for a, b in zip(vecs[i], vecs[j])]) for j in range(N)] for i in range(N)]
    smul = [[enc([(lam * a) % p for a in vecs[i]]) for i in range(N)] for lam in range(p)]
    unit = [enc([1 if k == i else 0 for k in range(n)]) for i in range(n)]
    # registers: 0..nx-1 x (N), nx..nx+ny-1 y (S), then slots, then copies
    fixed = {}
    for i in range(nx): fixed[i] = 0
    for i in range(nx, nx + ny): fixed[i] = 1
    flex = [i for i in range(R) if i not in fixed]
    start_c = tuple(unit[i] for i in range(n)) + tuple(0 for _ in range(ncopy))
    start_s = [0] * R
    for i in range(nx, nx + ny): start_s[i] = 1
    for k, (a, b) in enumerate(slots): start_s[nx + ny + k] = a
    # goal contents
    goal = {}
    for i in range(nx):
        if mode == 'GF': goal[i] = unit[i]
        else:
            v = [1 if k == i else 0 for k in range(n)]
            for T in range(ny):
                v[nx + T] = (v[nx + T] - P[T][i]) % p
            goal[i] = enc(v)
    for T in range(ny):
        if mode == 'GF':
            v = [1 if k == nx + T else 0 for k in range(n)]
            for i in range(nx): v[i] = (v[i] + P[T][i]) % p
            goal[nx + T] = enc(v)
        else: goal[nx + T] = unit[nx + T]
    for k in range(nq): goal[nx + ny + k] = unit[nx + ny + k]
    end_status = {nx + ny + k: b for k, (a, b) in enumerate(slots)}
    def bits(s): return sum(b << i for i, b in enumerate(s))
    s0 = (start_c, bits(start_s))
    # copies start at zero: their initial status is free -> start them in either status (zero content flips are free anyway)
    dist = {s0: 0}
    dq = deque([s0])
    best = None
    while dq:
        st = dq.popleft()
        c, sb = st
        dcur = dist[st]
        if best is not None and dcur >= best: continue
        if all(c[i] == g for i, g in goal.items()) and all(((sb >> i) & 1) == e for i, e in end_status.items()):
            best = dcur if best is None else min(best, dcur)
            continue
        # free moves
        for i in range(R):
            si = (sb >> i) & 1
            for j in range(R):
                if i == j or ((sb >> j) & 1) != si or c[j] == 0: continue
                for lam in range(1, p):
                    nc = list(c); nc[i] = add[c[i]][smul[lam][c[j]]]
                    ns = (tuple(nc), sb)
                    if ns not in dist or dist[ns] > dcur:
                        dist[ns] = dcur; dq.appendleft(ns)
        # status changes of flexible registers
        for i in flex:
            nsb = sb ^ (1 << i)
            w = 0 if c[i] == 0 else 1
            ns = (c, nsb)
            if ns not in dist or dist[ns] > dcur + w:
                dist[ns] = dcur + w
                if w: dq.append(ns)
                else: dq.appendleft(ns)
    return best, len(dist)


def rank_mod(P, p):
    M = [list(r) for r in P]; rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], p - 2, p)
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c] * inv % p
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


cases = [
    dict(name='A1 GF 1x1y, 1 copy', p=3, nx=1, ny=1, slots=[], copies=1, P=[[1]], mode='GF'),
    dict(name='A2 GB 1x1y, 1 copy', p=3, nx=1, ny=1, slots=[], copies=1, P=[[1]], mode='GB'),
    dict(name='A3 GF 1x1y, no copy, no slot', p=3, nx=1, ny=1, slots=[], copies=0, P=[[1]], mode='GF'),
    dict(name='B1 GF 1x1y, slot N->S (natural climb), no copy', p=3, nx=1, ny=1, slots=[(0, 1)], copies=0, P=[[1]], mode='GF'),
    dict(name='B2 GF 1x1y, slot N->S (natural climb), 1 copy', p=3, nx=1, ny=1, slots=[(0, 1)], copies=1, P=[[1]], mode='GF'),
    dict(name='B3 GF 1x1y, slot N->N, no copy', p=3, nx=1, ny=1, slots=[(0, 0)], copies=0, P=[[1]], mode='GF'),
    dict(name='B4 GF 1x1y, slot S->S, no copy', p=3, nx=1, ny=1, slots=[(1, 1)], copies=0, P=[[1]], mode='GF'),
    dict(name='B5 GB 1x1y, slot N->S, no copy', p=3, nx=1, ny=1, slots=[(0, 1)], copies=0, P=[[1]], mode='GB'),
    dict(name='B6 GB 1x1y, slot S->N, no copy', p=3, nx=1, ny=1, slots=[(1, 0)], copies=0, P=[[1]], mode='GB'),
    dict(name='C1 GF 2x2y P=I, 1 copy', p=2, nx=2, ny=2, slots=[], copies=1, P=[[1, 0], [0, 1]], mode='GF'),
    dict(name='C2 GF 2x2y P=[[1,1],[1,1]] (rank 1), 1 copy', p=2, nx=2, ny=2, slots=[], copies=1, P=[[1, 1], [1, 1]], mode='GF'),
    dict(name='C3 GB 2x2y P=I, 1 copy', p=2, nx=2, ny=2, slots=[], copies=1, P=[[1, 0], [0, 1]], mode='GB'),
    dict(name='C4 GF 2x2y P=[[1,1],[0,1]], 1 copy', p=2, nx=2, ny=2, slots=[], copies=1, P=[[1, 1], [0, 1]], mode='GF'),
]
out = []
for cs in cases:
    best, ns = run(cs)
    rk = rank_mod(cs['P'], cs['p'])
    flips = sum(1 for a, b in cs['slots'] if a != b)
    pred = rk + flips
    print('%-55s min paid status changes = %s | rank P + slot flips = %d | states %d' % (cs['name'], best, pred, ns))
    out.append(dict(case=cs['name'], min_changes=best, rank=rk, slot_flips=flips, lemma_lower=pred, states=ns))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
