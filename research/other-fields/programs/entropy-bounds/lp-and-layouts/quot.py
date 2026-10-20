# quot.py -- how much can a k-part contraction of Gamma_m certify at most?
# For a partition pi of the frames into k parts, the quotient network's routing value is
#   rho(pi) = E_{ordered opposite pairs}[d_q(pi(s), pi(t))] / (share of edges crossing pi)
# (d_q = distance in the quotient graph). rho(pi) is an UPPER limit on any coding bound the quotient can give
# (coding never costs more than routing); the cut LP is the max of rho over 2-part pi (= best cut).
# Also tests l1-embeddability of the quotient metric by LP (cut cone membership).
import sys, random, itertools, math, collections
import numpy as np
from scipy.optimize import linprog
import dpg

def quotient_metric(part, k, edges):
    nb = [set() for _ in range(k)]
    for (u, v) in edges:
        a, b = part[u], part[v]
        if a != b: nb[a].add(b); nb[b].add(a)
    INF = 10 ** 9
    D = [[INF] * k for _ in range(k)]
    for s in range(k):
        D[s][s] = 0; q = collections.deque([s])
        while q:
            x = q.popleft()
            for y in nb[x]:
                if D[s][y] == INF: D[s][y] = D[s][x] + 1; q.append(y)
    return D

def rho(part, k, edges, opp):
    D = quotient_metric(part, k, edges)
    cross = sum(1 for (u, v) in edges if part[u] != part[v])
    if cross == 0: return 0.0, D
    num = sum(D[part[s]][part[t]] for (s, t) in opp)
    if num >= 10 ** 9: return 0.0, D
    return (num / len(opp)) / (cross / len(edges)), D

def is_l1(D):
    k = len(D)
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    cuts = []
    for mask in range(1, 1 << (k - 1)):
        S = {i for i in range(k) if mask >> i & 1}
        cuts.append([1.0 if ((i in S) != (j in S)) else 0.0 for (i, j) in pairs])
    A = np.array(cuts).T; b = np.array([float(D[i][j]) for (i, j) in pairs])
    r = linprog(np.zeros(A.shape[1]), A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    return r.status == 0

def anneal(m, k, iters, seed, verts=None):
    rng = random.Random(seed)
    verts, idx, adj, dist = G[m]
    n = len(verts)
    edges = [(i, j) for i in range(n) for j in adj[i] if i < j]
    opp = [(i, j) for i in range(n) for j in range(n) if dist[i][j] == m]
    part = [rng.randrange(k) for _ in range(n)]
    cur, _ = rho(part, k, edges, opp); best = (cur, part[:])
    T = 0.05
    for it in range(iters):
        v = rng.randrange(n); old = part[v]; new = rng.randrange(k)
        if new == old: continue
        part[v] = new
        val, _ = rho(part, k, edges, opp)
        if val >= cur or rng.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0]: best = (val, part[:])
        else:
            part[v] = old
        T = max(0.002, T * 0.9995)
    return best, edges, opp

G = {}
if __name__ == '__main__':
    m = int(sys.argv[1]); ks = [int(x) for x in sys.argv[2].split(',')]; iters = int(sys.argv[3]); restarts = int(sys.argv[4])
    G[m] = dpg.build(m)
    for k in ks:
        res = []
        for r in range(restarts):
            (val, part), edges, opp = anneal(m, k, iters, 1000 * k + r)
            D = quotient_metric(part, k, edges)
            used = sorted(set(part))
            res.append((val, part, D))
        res.sort(key=lambda t: -t[0])
        val, part, D = res[0]
        sizes = collections.Counter(part)
        print(f"m={m} k={k}: best rho = {val:.5f} (top 3: {[round(r[0],5) for r in res[:3]]}); part sizes {sorted(sizes.values())}; "
              f"quotient metric l1? {is_l1([row[:] for row in D]) if len(D) <= 9 else 'n/a'}")
        print("   quotient distances:", D)
        sys.stdout.flush()
