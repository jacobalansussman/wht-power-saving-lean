"""cliquefloor.py (t-km): the counting floor of u-ours/toy.py (octads only: cst >= 53 on P112) redone with ALL maximal
cliques of the orthogonality graph on P112 (1260 octads, 512 7-cliques, 224 6-cliques).  Clique Hall needs every
maximal clique Q to meet >= |Q| distinct frames, so sum_j hits(F_j) >= sum_Q |Q| where hits(F) = #cliques meeting F.
H(d) = max hits over flats of dim d (exact, all flats from flats_cache.npz); floor = least sum d_j with
sum H(d_j) >= need (knapsack DP).  Also the weighted variant: weights w_Q free; we report octads-only (control 53),
all cliques, and 6/7-cliques only.  usage: python -I cliquefloor.py"""
import sys, itertools, json
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from model import *
import numpy as np
OM = [sum(1 << T for T in range(NP) if ORTH[S, T]) for S in range(NP)]
MC = []
def bk(Rr, Pm, Xm):
    if not Pm and not Xm: MC.append(Rr); return
    u = (Pm | Xm).bit_length() - 1; c = Pm & ~OM[u]
    while c:
        low = c & -c; v = low.bit_length() - 1; c ^= low
        bk(Rr + [v], Pm & OM[v], Xm & OM[v]); Pm &= ~low; Xm |= low
bk([], (1 << NP) - 1, 0)
def tomask(ports):
    a = np.zeros(NP, bool); a[list(ports)] = True
    return np.packbits(np.concatenate([a, np.zeros(128 - NP, bool)])).view('>u8').astype(np.uint64)
CQ = np.array([tomask(q) for q in MC]); size = np.array([len(q) for q in MC])
cz = np.load('flats_cache.npz')
H = {}
for k in range(1, 9):
    P8 = cz['d%d' % k]
    P16 = np.concatenate([P8, np.zeros((len(P8), 2), np.uint8)], axis=1).view('>u8').astype(np.uint64)  # n x 2
    best = {'all': 0, 'oct': 0, 'small': 0}
    for a in range(0, len(P16), 20000):
        X = P16[a:a + 20000]
        meet = ((X[:, None, 0] & CQ[None, :, 0]) | (X[:, None, 1] & CQ[None, :, 1])) != 0     # chunk x 1996
        best['all'] = max(best['all'], int(meet.sum(1).max()))
        best['oct'] = max(best['oct'], int(meet[:, size == 8].sum(1).max()))
        best['small'] = max(best['small'], int(meet[:, size < 8].sum(1).max()))
        # weighted by |Q| (same need side): hits weighted
    H[k] = best
    print('dim', k, 'flats', len(P8), 'max cliques met', best, flush=True)
def floor(key, need):
    dp = [0] + [-1] * 300
    for c in range(1, 301):
        dp[c] = max(dp[c - k] + H[k][key] for k in H if k <= c)
    return next(c for c in range(301) if dp[c] >= need)
need = {'oct': int((size == 8).sum()) * 8, 'small': int(size[size < 8].sum()), 'all': int(size.sum())}
res = {key: dict(need=need[key], floor=floor(key, need[key])) for key in need}
print('counting floors (cst >=):', res)
json.dump(dict(H={str(k): v for k, v in H.items()}, floors=res, cliques={str(s): int((size == s).sum()) for s in (6, 7, 8)}),
          open('cliquefloor.json', 'w'), indent=1)
