#!/usr/bin/env python3
"""lagfloor.py (t-family) -- the exact copy-floor toy of ../lib/copyfloor.py, extended to the
ENLARGED (twisted, directed) knows relation of Lagrangian label families.  Standard library only; python3 -I.

A family with fixed twists: port T has source vector p_T and target vector q_T (pulled back into one symplectic
space F_2^{2h}), omega(p_T, q_T) = 1.  N[T][S] = omega(p_T, q_S).  Target S knows source T iff N[T][S] = 0
(necessary condition; the exact frame condition is checked separately in lagframes.py).  The hub's labels are the
case p_T = (u_T|u_T), q_T = (0|u_T): N = u_T.u_S, symmetric.
A decoder over the payload field F_p is a list of terms (r_j, g_j); B = sum g_j r_j^T must have B[S][S] != 0 and
B[S][T] = 0 whenever S does not know T.  A term reading R and serving D costs (twisted pairing floor, t-sandwich)
pair(R, D) = rank_F2 N[R, D].
Field-match lemma [proved on paper]: if p = 2 then sum_j pair_j >= rank_F2 of B o N^T = diag(B) = v, so x >= 1 for
every family and every twist; a saving needs the payload characteristic (0 on the complex side) to differ from the
frame field F_2.  The toy therefore uses p = 3 and p = 5 as stand-ins for characteristic 0.

  python3 -I lagfloor.py all4 [p]     every unit-diagonal F_2 pattern N on v = 4 ports (up to relabelling): floors
  python3 -I lagfloor.py rand v n [p] n random patterns on v ports
"""
import heapq, itertools, sys, os, random, time

def rank2(rows):
    basis = []
    for x in rows:
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    return len(basis)

class DUnit:
    """directed unit from a 0/1 matrix N (list of lists, N[T][S]); costs = F_2 rank of N[R, D]."""
    def __init__(self, N, p=3, cost=None):
        self.N, self.v, self.p = N, len(N), p
        assert all(N[i][i] == 1 for i in range(self.v))
        self.knows = [[S != T and N[T][S] == 0 for T in range(self.v)] for S in range(self.v)]
        self._c = {}
        self.costfn = cost
    def pair(self, R, D):
        key = (R, D)
        if key not in self._c:
            if self.costfn is not None:
                self._c[key] = self.costfn(R, D)
            else:
                rows = [sum(self.N[T][S] << i for i, S in enumerate(D)) for T in R]
                self._c[key] = rank2(rows)
        return self._c[key]

def exact_floor(U, limit_states=3_000_000):
    """Dijkstra over the constrained entries of B (as copyfloor.exact_floor), cost = U.pair(R, D)."""
    p, v = U.p, U.v
    E = [(S, S) for S in range(v)] + [(S, T) for S in range(v) for T in range(v) if S != T and not U.knows[S][T]]
    idx = {e: i for i, e in enumerate(E)}; nE = len(E)
    def vecs(normalise):
        out = []
        for t in itertools.product(range(p), repeat=v):
            if any(t):
                if normalise and t[next(i for i in range(v) if t[i])] != 1:
                    continue
                out.append(t)
        return out
    best = {}
    for r in vecs(True):
        Rs = tuple(T for T in range(v) if r[T])
        for g in vecs(False):
            Ds = tuple(S for S in range(v) if g[S])
            eff = [0] * nE; nz = False
            for S in Ds:
                for T in Rs:
                    i = idx.get((S, T))
                    if i is not None:
                        eff[i] = g[S] * r[T] % p; nz = nz or eff[i] != 0
            if not nz:
                continue
            c = U.pair(Rs, Ds)
            if c is None:
                continue
            key = tuple(eff)
            if key not in best or c < best[key][0]:
                best[key] = (c, (r, g))
    moves = sorted(((c, eff, tg) for eff, (c, tg) in best.items()), key=lambda m: m[0])
    def goal(st):
        return all(st[i] % p for i in range(v)) and all(st[i] % p == 0 for i in range(v, nE))
    start = tuple([0] * nE); dist = {start: 0}; prev = {}; pq = [(0, start)]; n = 0
    while pq:
        d, st = heapq.heappop(pq)
        if d > dist.get(st, 1 << 30):
            continue
        if goal(st):
            path = []
            while st in prev:
                st, tg = prev[st]; path.append(tg)
            return d, path[::-1], n
        n += 1
        if n > limit_states:
            return None, None, n
        for c, eff, tg in moves:
            ns = tuple((a + b) % p for a, b in zip(st, eff)); nd = d + c
            if nd < dist.get(ns, 1 << 30):
                dist[ns] = nd; prev[ns] = (st, tg); heapq.heappush(pq, (nd, ns))
    return None, None, n

def canon(N):
    v = len(N)
    best = None
    for pm in itertools.permutations(range(v)):
        key = tuple(N[pm[i]][pm[j]] for i in range(v) for j in range(v))
        if best is None or key < best:
            best = key
    return best

def all_patterns(v):
    off = [(i, j) for i in range(v) for j in range(v) if i != j]
    seen = set()
    for bits in range(1 << len(off)):
        N = [[int(i == j) for j in range(v)] for i in range(v)]
        for k, (i, j) in enumerate(off):
            if bits >> k & 1:
                N[i][j] = 1
        c = canon(N)
        if c in seen:
            continue
        seen.add(c)
        yield N

if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all4'
    HERE = os.path.dirname(os.path.abspath(__file__))
    if cmd == 'all4':
        p = int(sys.argv[2]) if len(sys.argv) > 2 else 3
        v = 4; hist = {}; best = None; t = time.time(); cnt = 0
        for N in all_patterns(v):
            U = DUnit(N, p)
            f, path, n = exact_floor(U)
            cnt += 1
            hist[f] = hist.get(f, 0) + 1
            if best is None or f < best[0]:
                best = (f, N, path)
        print('v=4, payload F_%d: %d patterns up to relabelling; floor histogram %s; min %s (%.1fs)'
              % (p, cnt, dict(sorted(hist.items(), key=lambda kv: (kv[0] is None, kv[0] or 0))), best[0], time.time() - t))
        print('  a minimiser N =', best[1], 'decoder', best[2])
    elif cmd == 'rand':
        v, n = int(sys.argv[2]), int(sys.argv[3]); p = int(sys.argv[4]) if len(sys.argv) > 4 else 3
        rng = random.Random(1); hist = {}; best = None; t = time.time()
        for _ in range(n):
            N = [[1 if i == j else rng.randint(0, 1) for j in range(v)] for i in range(v)]
            f, path, nn = exact_floor(DUnit(N, p))
            hist[f] = hist.get(f, 0) + 1
            if f is not None and (best is None or f < best[0]):
                best = (f, N, path)
        print('v=%d, payload F_%d: %d random patterns; floor histogram %s (%.1fs)' % (v, p, n, hist, time.time() - t))
        print('  a minimiser N =', best[1], 'decoder', best[2])
