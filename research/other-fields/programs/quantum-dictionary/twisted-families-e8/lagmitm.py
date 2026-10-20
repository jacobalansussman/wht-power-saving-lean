#!/usr/bin/env python3
"""lagmitm.py (t-family) -- exact test "floor <= v - 1 ?" for small directed units, by meet in the middle.
Same model as lagfloor.py (term cost = rank_F2 N[R, D], payload F_p).  A decoder of cost <= v-1 exists iff some
x in S_a and y in S_b (S_c = effects on the constrained entries reachable with total cost <= c) with a + b <= v-1
have offdiag(x) + offdiag(y) = 0 and every diagonal entry of x + y non-zero.  Floor = v otherwise (one term per
port always works).  Standard library only.
usage: python3 -I lagmitm.py all v p      every unit-diagonal pattern on v ports up to relabelling (v <= 4)
       python3 -I lagmitm.py rand v n p   n random patterns on v ports (v = 5)"""
import itertools, sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lagfloor import rank2, all_patterns

def terms(N, p):
    v = len(N)
    knows = [[S != T and N[T][S] == 0 for T in range(v)] for S in range(v)]
    off = [(S, T) for S in range(v) for T in range(v) if S != T and not knows[S][T]]
    oi = {e: i for i, e in enumerate(off)}
    def vecs(norm):
        for t in itertools.product(range(p), repeat=v):
            if any(t) and (not norm or t[next(i for i in range(v) if t[i])] == 1):
                yield t
    G = list(vecs(False)); Rv = list(vecs(True))
    best = {}
    for r in Rv:
        Rs = [T for T in range(v) if r[T]]
        for g in G:
            Ds = [S for S in range(v) if g[S]]
            c = rank2([sum(N[T][S] << i for i, S in enumerate(Ds)) for T in Rs])
            od = [0] * len(off)
            for S in Ds:
                for T in Rs:
                    i = oi.get((S, T))
                    if i is not None: od[i] = g[S] * r[T] % p
            dg = tuple(g[S] * r[S] % p for S in range(v))
            key = (tuple(od), dg)
            if key not in best or c < best[key]:
                best[key] = c
    return best, len(off)

def floor_lt_v(N, p):
    v = len(N)
    best, noff = terms(N, p)
    budget = v - 1
    # S[c] : dict effect -> True for min cost exactly <= c  (effect = (offdiag tuple, diag tuple))
    lay = {(tuple([0] * noff), tuple([0] * v)): 0}
    by_cost = {}
    for k, c in best.items():
        by_cost.setdefault(c, []).append(k)
    half = budget // 2
    frontier = dict(lay)
    reach = dict(lay)
    # expand up to cost ceil(budget/2) by adding terms
    hi = budget - half
    changed = True
    while changed:
        changed = False
        new = {}
        for e, ce in reach.items():
            for c, ks in by_cost.items():
                if ce + c > hi: continue
                for k in ks:
                    od = tuple((a + b) % p for a, b in zip(e[0], k[0])); dg = tuple((a + b) % p for a, b in zip(e[1], k[1]))
                    ne = (od, dg)
                    if reach.get(ne, 1 << 30) > ce + c and new.get(ne, 1 << 30) > ce + c:
                        new[ne] = ce + c
        for ne, c in new.items():
            if reach.get(ne, 1 << 30) > c:
                reach[ne] = c; changed = True
    idx = {}
    for e, c in reach.items():
        idx.setdefault(e[0], []).append((c, e[1]))
    for e, c in reach.items():
        if c > half: continue
        neg = tuple((-a) % p for a in e[0])
        for c2, dg2 in idx.get(neg, ()):
            if c + c2 <= budget and all((a + b) % p for a, b in zip(e[1], dg2)):
                return True, c + c2
    return False, v

if __name__ == '__main__':
    cmd, v = sys.argv[1], int(sys.argv[2])
    if cmd == 'all':
        p = int(sys.argv[3]); t = time.time(); cnt = hits = 0; ex = None
        for N in all_patterns(v):
            ok, c = floor_lt_v(N, p); cnt += 1
            if ok:
                hits += 1; ex = ex or (N, c)
        print('v=%d payload F_%d: %d patterns up to relabelling, %d with floor < v %s (%.1fs)'
              % (v, p, cnt, hits, ('e.g. N=%s cost %d' % ex) if ex else '', time.time() - t))
    elif cmd == 'rand':
        n, p = int(sys.argv[3]), int(sys.argv[4]); rng = random.Random(7); t = time.time(); hits = 0; ex = None
        for _ in range(n):
            N = [[1 if i == j else rng.randint(0, 1) for j in range(v)] for i in range(v)]
            ok, c = floor_lt_v(N, p)
            if ok:
                hits += 1; ex = ex or (N, c)
        print('v=%d payload F_%d: %d random patterns, %d with floor < v %s (%.1fs)'
              % (v, p, n, hits, ('e.g. N=%s cost %d' % ex) if ex else '', time.time() - t))
