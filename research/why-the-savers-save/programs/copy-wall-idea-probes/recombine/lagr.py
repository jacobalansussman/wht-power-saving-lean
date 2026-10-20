"""lagr.py (g-recombine) -- do Lagrangian (Clifford) frames give a cheaper single-total copy than subspace frames?
Frames = Lagrangians of F_2^{2m} (x|z), symplectic form x.z' + x'.z; subspace frame F <-> L_F = F (+) F^perp.
d(L,L') = m - dim(L cap L').  (1) interval check: every Lagrangian on a geodesic L_A..L_B (A <= B) is a subspace
frame L_C with A <= C <= B.  (2) Steiner check: one total reading sources with labels a_1..a_k; source T is free
anywhere on ITS climb chain (we allow any chain: every subspace frame containing a_T); the total must reach a
target cut frame L_0 (frame 0).  Min Steiner tree in the dual polar graph vs dim span(a).  Exhaustive, m = 3.
usage: python3 -I lagr.py [m]"""
import itertools, sys
m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
n = 2 * m
def omega(u, w):
    x1, z1 = u >> m, u & ((1 << m) - 1); x2, z2 = w >> m, w & ((1 << m) - 1)
    return (bin(x1 & z2).count('1') + bin(x2 & z1).count('1')) & 1
def span(vs):
    S = {0}
    for v in vs:
        S |= {s ^ v for s in S}
    return frozenset(S)
# enumerate Lagrangians: isotropic m-dim subspaces
lag = set()
def rec(basis, S, start):
    if len(basis) == m:
        lag.add(S); return
    for v in range(start, 1 << n):
        if v in S: continue
        if all(omega(v, b) == 0 for b in basis) and omega(v, v) == 0:
            rec(basis + [v], frozenset(S | {s ^ v for s in S}), v + 1)
rec([], frozenset({0}), 1)
lag = list(lag); idx = {L: i for i, L in enumerate(lag)}
print('m', m, 'Lagrangians', len(lag))
def dim(S): return len(S).bit_length() - 1
def dist(L, M): return m - dim(L & M)
N = len(lag)
D = [[dist(lag[i], lag[j]) for j in range(N)] for i in range(N)]
mask = (1 << m) - 1
def LF(F):  # F = frozenset of x-vectors (subspace of F_2^m); L_F = {(x, z): x in F, z in F^perp}
    perp = [z for z in range(1 << m) if all(bin(z & x).count('1') % 2 == 0 for x in F)]
    return frozenset((x << m) | z for x in F for z in perp)
subs = set()
for k in range(m + 1):
    for vs in itertools.combinations(range(1, 1 << m), k):
        subs.add(span(vs))
subs = list(subs)
subframe = {idx[LF(F)]: F for F in subs}
print('subspace frames', len(subframe), 'non-subspace Lagrangians', N - len(subframe))
# (1) interval check
bad = 0
for A in subs:
    for B in subs:
        if A <= B:
            ia, ib = idx[LF(A)], idx[LF(B)]
            for c in range(N):
                if D[ia][c] + D[c][ib] == D[ia][ib]:
                    if c not in subframe or not (A <= subframe[c] <= B):
                        bad += 1
print('(1) interval check: Lagrangians on a geodesic L_A..L_B (A<=B) that are not subspace frames between:', bad)
# (2) Steiner: terminals = groups (source T: all subspace frames containing a_T) + target L_0.
# Group Steiner tree via Dreyfus-Wagner over groups (k <= 4).
import heapq
def group_steiner(groups):
    k = len(groups); full = (1 << k) - 1
    INF = 10 ** 9
    dp = [[INF] * N for _ in range(1 << k)]
    for g, G in enumerate(groups):
        for v in G: dp[1 << g][v] = 0
    for S in range(1, full + 1):
        sub = (S - 1) & S
        while sub:
            if sub < S ^ sub:
                for v in range(N):
                    c = dp[sub][v] + dp[S ^ sub][v]
                    if c < dp[S][v]: dp[S][v] = c
            sub = (sub - 1) & S
        # relax by shortest paths (unit edges: distance-1 neighbours)
        best = dp[S]
        for v in range(N):
            for u in range(N):
                c = best[u] + D[u][v]
                if c < best[v]: best[v] = c
    return min(dp[full])
L0 = idx[LF(frozenset({0}))]
worst_gap = []
for k in (2, 3, 4):
    for labs in itertools.combinations(range(1, 1 << m), k):
        if k > m + 1: break
        groups = [[i for i, F in subframe.items() if a in F] for a in labs] + [[L0]]
        st = group_steiner(groups)
        dP = dim(span(labs))
        worst_gap.append((st - dP, labs, st, dP))
print('(2) group-Steiner minus dim P over label sets (k = 2..4): min gap', min(g[0] for g in worst_gap),
      ' max gap', max(g[0] for g in worst_gap), ' cases', len(worst_gap))
print('    examples', sorted(worst_gap)[:3])
