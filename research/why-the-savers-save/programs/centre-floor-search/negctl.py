#!/usr/bin/env python3
"""ctr-search / negctl.py: negative controls for the checker.  (1) the hub's 64 decoder with every gather support cut
to a 6-dim label subspace (cost 48 < 52.36): must fail the exact check; (2) a fabricated 'decoder' that is exact
but whose cost is under-reported would be REFUSED by dsearch.verify's floor guard: simulated by calling the guard
on the frame decoder with the floor raised to 70 (it must refuse)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p112 import *
import dsearch
from calib import frame_decoder
F1 = [tuple(1 if j == i else 0 for j in range(8)) for i in range(8)]
Phi, Psi = frame_decoder(F1)
# (1) cut each gather support to the ports whose labels lie in a 6-dim subspace of its 8-dim span
from nd import rref_subspaces, mask
cut = []
for k in range(8):
    A = [b for b in range(N) if Phi[k][b] != 0]
    best = 0; bm = 0
    for rows in rref_subspaces(9, 3):
        m = (1 << N) - 1
        for r in rows: m &= mask[r]
        c = sum(1 for b in A if m >> b & 1)
        if c > best: best, bm = c, m
    cut.append([Phi[k][b] if (bm >> b & 1) else 0 for b in range(N)])
ok, bad = check_exact(cut, Psi); c = costs(cut, Psi)
print('(1) frame decoder cut to 6-dim gathers: one-sided cost', c['star_gather'], 'exact', ok, 'violations', bad, '-> refused' if not ok else '-> ACCEPTED (bug!)')
# (2) floor guard
used = list(range(8)); reps = [[tuple(i for i in range(8) if LINES[b][i] != 0)] for b in range(N)]
dsearch.FLOOR_ONE_SIDED = 70.0
v = dsearch.verify(F1, reps, [0] * N, 'frame with floor 70')
print('(2) guard with floor 70 on the 64 decoder:', v.get('REFUSED', 'not refused (bug!)'))
