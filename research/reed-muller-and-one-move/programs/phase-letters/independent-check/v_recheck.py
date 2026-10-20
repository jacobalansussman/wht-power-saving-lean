"""Re-test rm2-phase's survivor lists with MY D + C + Z (v_cls leaf logic) and count tight forms.
Targets are read through rm2-phase's own index -> matrix map (tile_alt.target_matrices), directions as ints (bit i = coord i)."""
import sys, numpy as np
sys.path.insert(0, '../ends-only-block4')
from tile_alt import target_matrices
from vb import *
mats, nonalt = target_matrices(4)
m = 4; Qz = Qz_matrix(m); nf, nd = Qz.shape
def maxclique(cand, adj):
    best = 0
    def bk(R, P, X):
        nonlocal best
        if P == 0 and X == 0:
            best = max(best, R); return
        if R + bin(P).count('1') <= best: return
        u = (P | X) & -(P | X); ui = u.bit_length() - 1
        Q = P & ~adj[ui]
        while Q:
            v = Q & -Q; vi = v.bit_length() - 1
            bk(R + 1, P & adj[vi], X & adj[vi]); P &= ~v; X |= v; Q &= ~v
    bk(0, cand, 0)
    return best
def test(dirs, tg):
    W = len(tg); inst = [z - 1 for z in dirs]; n = len(inst)
    cnt = np.zeros(nd, dtype=np.int64)
    for d in inst: cnt[d] += 1
    Dq = cnt @ Qz.T
    r = pair_matrix(m, [np.array(mats[t]) for t in tg]).sum(1)
    if (Dq < r).any(): return 'D', None
    fm = [sum(1 << k for k, d in enumerate(inst) if Qz[f, d]) for f in range(nf)]
    tightf = np.where((Dq == r) & (r >= 1))[0]
    adj = [0] * n; live = 0
    for f in tightf:
        mk = fm[f]; live |= mk
        for k in range(n):
            if mk >> k & 1: adj[k] |= mk & ~(1 << k)
    c = maxclique(live, adj)
    if c > W: return 'C', len(tightf)
    for f in np.where(Dq + r < 2 * c)[0]:
        if Dq[f] + r[f] < 2 * maxclique(fm[f] & live, adj): return 'Z', len(tightf)
    return 'pass', len(tightf)
import re
res = {}
for line in open('../ends-only-block4/pass414.txt'):
    mm = re.search(r"dirs \(([\d, ]+)\) targets \(([\d, ]+)\)", line)
    d = [int(x) for x in mm.group(1).split(',')]; t = [int(x) for x in mm.group(2).split(',')]
    v, nt = test(d, t); res.setdefault(('414', v, nt > 0, sum(x >= nonalt for x in t) > 0), []).append(t)
step = int(sys.argv[1]) if len(sys.argv) > 1 else 1
for i, line in enumerate(open('../ends-only-block4/surv_fe_alt_415_all.txt')):
    if i % step: continue
    a, b = line.strip().split('|')
    d = [int(x) for x in a.split(',')]; t = [int(x) for x in b.split(',')]
    v, nt = test(d, t); res.setdefault(('415', v, nt > 0, sum(x >= nonalt for x in t)), []).append(t)
for k in sorted(res, key=str): print(k, len(res[k]))
