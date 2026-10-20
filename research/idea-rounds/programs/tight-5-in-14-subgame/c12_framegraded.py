# c12_framegraded.py -- "frame-graded hooks": Theorem R on V56 (x) C[G/B] (diagonal action) with the B-invariant
# certificate Y' = sum_P Z_{d(P)} (x) e_P, Z_d in {0, Y, Yperp, V56} chosen per frame distance d = dist(L_Z, P).
# cap(h)  = sum_P [dim(Z_{d(P)} + h Z_{d(h^-1 P)}) - dim Z_{d(P)}]
# dem(K)  = sum_P [dim(Z_{d(P)} + K Z_{d(K^-1 P)}) - dim Z_{d(P)}]
# Calibration: all Z_d = Y gives 1350 / 3780 = 2.8 (Lemma A); question: does any choice beat 2.8?
# Also the subgame version (30 singular frames, free group B' keeps every graph-certificate Y(t), and Z_d may vary in
# the family per distance).
import itertools, numpy as np
from collections import Counter
from lib514 import *
R = Irreps()
dat = np.load('v56_basis.npz'); Vw = dat['V']; Y = Vw[:, dat['o28']]; Yp = Vw[:, dat['o7']]
lag = set()
for vs in itertools.combinations(range(1, NV), 3):
    S = span(vs)
    if len(S) == 8 and all(BTAB[a, b] == 0 for a in S for b in S): lag.add(S)
lag = sorted(lag, key=sorted)
def dist(F1, F2): return 3 - {1: 0, 2: 1, 4: 2, 8: 3}[len(F1 & F2)]
def img(g, F): return frozenset(int(g[v]) for v in F)
hinv = inv(h1); Kinv = inv(K)
nh = Counter((dist(LZ, P), dist(LZ, img(hinv, P))) for P in lag)
nK = Counter((dist(LZ, P), dist(LZ, img(Kinv, P))) for P in lag)
print('h: (d(P), d(h^-1 P)) counts', dict(nh)); print('K: (d(P), d(K^-1 P)) counts', dict(nK))
subs = {'0': np.zeros((56, 0)), 'Y': Y, 'P': Yp, 'V': np.eye(56)}
def dsum(Z, g, Z2):
    if Z.shape[1] + Z2.shape[1] == 0: return 0
    return rank(np.hstack([Z, R.r56(g) @ Z2]))
tab = {}
for a, b in itertools.product(subs, repeat=2):
    tab[(a, b, 'h')] = dsum(subs[a], h1, subs[b]) - subs[a].shape[1]
    tab[(a, b, 'K')] = dsum(subs[a], K, subs[b]) - subs[a].shape[1]
best = []
for choice in itertools.product(subs, repeat=4):
    cap = sum(n * tab[(choice[d], choice[e], 'h')] for (d, e), n in nh.items())
    dem = sum(n * tab[(choice[d], choice[e], 'K')] for (d, e), n in nK.items())
    if cap > 0: best.append((dem / cap, cap, dem, ''.join(choice)))
best.sort(reverse=True)
print('best frame-graded choices (Z_0 Z_1 Z_2 Z_3 with 0/Y/P=Yperp/V):')
for b in best[:8]: print('   ratio %.4f cap %d demand %d  %s' % b)
print('all-Y calibration:', [b for b in best if b[3] == 'YYYY'])
