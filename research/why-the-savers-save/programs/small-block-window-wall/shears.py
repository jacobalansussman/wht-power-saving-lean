#!/usr/bin/env python3
"""an2-smallnet / shears.py (own code, stdlib).  Question: in a crossing network built from bank-to-bank shear calls
(b_i += c * b_j, c = +-1, the hub's GF / GB calls), every data bank must end holding another bank's data (a
fixed-point-free signed permutation of banks; otherwise the bank gains nothing).  Each call a bank takes part in costs
it one window of k-1 fresh frame levels (the stage-window wall), so the block is at least 1 + (k-1) * S where S is the
largest number of calls any bank takes part in.  This program finds, for n = 2..6 banks, the least S and the least
number of calls, by breadth-first search over shear sequences (integer matrices, exact)."""
import itertools, json, os, sys
from collections import deque

def is_signed_derangement(M, n):
    for i in range(n):
        row = M[i]
        nz = [j for j in range(n) if row[j] != 0]
        if len(nz) != 1 or abs(row[nz[0]]) != 1 or nz[0] == i:
            return False
    cols = sorted(next(j for j in range(n) if M[i][j]) for i in range(n))
    return cols == list(range(n))

def search(n, maxlen, maxper):
    """BFS over states (matrix, per-bank call counts); returns the shortest sequence reaching a signed derangement
    with every bank in at most maxper calls, or None."""
    I = tuple(tuple(1 if i == j else 0 for j in range(n)) for i in range(n))
    start = (I, (0,) * n)
    seen = {start}
    q = deque([(start, ())])
    while q:
        (M, cnt), seq = q.popleft()
        if len(seq) >= maxlen:
            continue
        for i in range(n):
            for j in range(n):
                if i == j or cnt[i] >= maxper or cnt[j] >= maxper:
                    continue
                for c in (1, -1):
                    R = [list(r) for r in M]
                    R[i] = [R[i][k] + c * M[j][k] for k in range(n)]
                    if max(abs(v) for r in R for v in r) > 3:
                        continue
                    N = tuple(tuple(r) for r in R)
                    nc = list(cnt); nc[i] += 1; nc[j] += 1
                    st = (N, tuple(nc))
                    if st in seen:
                        continue
                    seen.add(st)
                    s2 = seq + ((i, j, c),)
                    if is_signed_derangement(N, n):
                        return s2
                    q.append((st, s2))
    return None

out = {}
for n in range(2, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 6):
    res = {}
    for maxper in (2, 3):
        maxlen = n * maxper // 2
        s = search(n, maxlen, maxper)
        res[maxper] = s
        print(f"n={n} banks, every bank in <= {maxper} calls: {'POSSIBLE with ' + str(len(s)) + ' calls ' + str(s) if s else 'impossible'}", flush=True)
    out[n] = {str(k): v for k, v in res.items()}
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'shears.json'), 'w'), indent=1)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'shears.json'), 'w'), indent=1)
