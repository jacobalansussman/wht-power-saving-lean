#!/usr/bin/env python3
"""an-ours: the per-replica rank identity of the P112 nine-stage network (own arithmetic; inputs: fsp-price/data/words.json
per-call profiles and widths, x-hosting/tools/host_price.py's child table CH, PAIR (our own tools, read as text))."""
import math, json
M = 73
words = {
 '588atm3': ({1:2484, 2:492, 3:280, 4:232, 5:56, 6:56, 8:8}, {4:48, 5:8, 7:32, 8:140, 9:276}),
 'U8E1': ({1:2548, 2:338, 3:446, 4:162, 5:72, 8:8}, {4:80, 5:32, 7:54, 8:122, 9:216}),
 'U8XM_5cc': ({1:2553, 2:338, 3:446, 4:162, 5:72, 8:13}, {4:80, 5:32, 7:54, 8:122, 9:221}),
 'T74 gauge44 (from T74E = 72 local + bank rest)': ({1:2340, 2:508, 3:420, 4:126, 5:82, 8:8}, None),
}
rk = lambda H: sum(r * n for r, n in H.items())
nm = lambda H: sum(r * n * math.log(M / r) for r, n in H.items())
CH = [('A2', 16, 2), ('B2', 8, 2), ('A3', 32, 3), ('B3', 8, 3), ('A4', 48, 4), ('B4', 8, 4), ('F1', 48, 1), ('F2', 32, 2), ('F3', 16, 3)]
PAIR = [1, 1, 1, 2, 2, 3, 3, 4, 4]
V, SIDES, STAGE_STEP = 112, 2, 8
print('per data array of pair p: in-call levels (8 per active stage) + idle children = paid levels')
for p in (1, 2, 3, 4):
    incall = STAGE_STEP * PAIR.count(p); idle = sum(r for _, r, q in CH if q == p)
    print(f'  pair {p}: {incall} in-call + {idle} idle = {incall + idle} of {M}  (free: {M - incall - idle})')
children = SIDES * V * sum(r for _, r, _ in CH)
data_arrays = SIDES * V * 4
out = {}
for k, (loc, wid) in words.items():
    L = rk(loc); N = nm(loc)
    line = f'{k:48s} per call: rank {L}, narrow {N:9.1f}, rank-1 share of narrow {loc[1]*math.log(M)/N:.3f}'
    if wid:
        Rw = rk(wid)
        per_rep = 9 * L + children - 9 * Rw
        D = M * data_arrays - per_rep
        line += f'; residual widths {Rw}; local - residuals = {L - Rw} (= 224*8 + 64 = {224*8+64}); D per replica = {D}'
        out[k] = dict(rank=L, narrow=N, residuals=Rw, local_minus_res=L - Rw, D_per_replica=D)
    else:
        out[k] = dict(rank=L, narrow=N)
    print(line)
print('data arrays per replica', data_arrays, '; children rank per replica', children, '; 73*896 - 320 =', 73 * 896 - 320)
json.dump(out, open('identity.json', 'w'), indent=1)
