"""small_play.py (g-small): no-hypothesis play on tiny units.
 (1) random bit-type units (Q-lines in {-1,0,1}^h), v <= 5: star vs pairing floor, payload F_2 vs F_3;
 (2) VECTOR-LINEAR with block-invertible diagonal blocks (t symbols per port, same label; the t x t block B_SS only
     has to be invertible, as a target can undo a local t x t transform on its own symbols) vs copyfloor's
     diagonal-dup model.  Report every case where a floor is below v (or t v), or pair < star, or p3 < p2."""
import sys, itertools, random, heapq
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from common import cf

def floor_blockinv(U, t, kind="star", limit=150000):
    """Dijkstra like cf.exact_floor, payload F_2, ports = U.v // t groups of t consecutive symbols with one label."""
    v = U.v; n = v // t; port = [i // t for i in range(v)]
    E = []
    for P in range(n):
        for a in range(t):
            for b in range(t): E.append((P * t + a, P * t + b))
    for S in range(v):
        for T in range(v):
            if port[S] != port[T] and not U.knows[S][T]: E.append((S, T))
    idx = {e: i for i, e in enumerate(E)}; nE = len(E); nb = n * t * t
    best = {}
    for r in range(1, 1 << v):
        Rs = [T for T in range(v) if r >> T & 1]
        for g in range(1, 1 << v):
            Ds = [S for S in range(v) if g >> S & 1]
            eff = 0
            for S in Ds:
                for T in Rs:
                    i = idx.get((S, T))
                    if i is not None: eff |= 1 << i
            if not eff: continue
            c = U.star(Rs) if kind == 'star' else U.pair(Rs, Ds)
            if eff not in best or c < best[eff]: best[eff] = c
    moves = sorted((c, e) for e, c in best.items())
    def goal(st):
        if st >> nb: return False
        for P in range(n):
            M = [[st >> (P * t * t + a * t + b) & 1 for b in range(t)] for a in range(t)]
            if cf.rank(M, 2) < t: return False
        return True
    dist = {0: 0}; pq = [(0, 0)]; k = 0
    while pq:
        d, st = heapq.heappop(pq)
        if d > dist[st]: continue
        if goal(st): return d, k
        k += 1
        if k > limit: return None, k
        for c, e in moves:
            ns = st ^ e; nd = d + c
            if nd < dist.get(ns, 1 << 30):
                dist[ns] = nd; heapq.heappush(pq, (nd, ns))
    return None, k

def main():
    out = []
    rng = random.Random(7)
    vecs3 = [list(x) for x in itertools.product((-1, 0, 1), repeat=3) if any(x)]
    lines3 = []
    for x in vecs3:
        if [-y for y in x] not in lines3: lines3.append(x)
    anomalies = 0; tried = 0
    for trial in range(12):
        v = 4
        labels = rng.sample(lines3, v)
        U2 = cf.Unit('r', 'Q', labels, p=2)
        if sum(map(sum, U2.knows)) == 0: continue
        tried += 1
        s2 = cf.star_floor_f2(U2)[0]
        p2 = cf.exact_floor(U2, "pair", limit_states=150000)[0]
        U3 = cf.Unit('r', 'Q', labels, p=3)
        s3 = None
        line = 'v %d knows %2d  star2 %s pair2 %s star3 %s  labels %s' % (v, sum(map(sum, U2.knows)), s2, p2, s3, labels)
        if s2 != v or p2 != v or (s3 is not None and s3 != v):
            anomalies += 1; line = 'ANOMALY ' + line
        out.append(line); print(line, flush=True)
    out.append('random bit-type units tried %d, anomalies %d' % (tried, anomalies))
    # vector-linear, block-invertible, t = 2
    for name in ['A2']:
        U = cf.make(name, 2, 2)
        if U.v > 8: continue
        for kind in ('star', 'pair'):
            d, k = floor_blockinv(U, 2, kind)
            print('block-invertible t=2 %-9s %s floor %s (2v = %d, states %d)' % (name, kind, d, U.v, k), flush=True)
    for o in out: print(o)
main()
