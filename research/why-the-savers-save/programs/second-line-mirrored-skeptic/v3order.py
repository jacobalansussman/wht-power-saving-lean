#!/usr/bin/env python3
"""v-cx-mirror / v3order.py (own code, stdlib).  The claim's own first-order cleaning model, re-searched.

Model (exactly cx-mirror's chaincost model): one register per total; it meets x_b only when its frame F lies in
v_b^perp (equivalently, F (inside the window) is orthogonal to g_b); meeting adds u_b to F; climbs are free; a drop
F -> F' costs dim F - dim F'.  chaincost.py drops, before each maximum run, to the meet of the WHOLE next run's
hyperplanes (cost up to |run|), with greedy run choice.  Here I drop LAZILY, one x at a time: if F is not inside
v_b^perp, drop to F & v_b^perp, which always costs exactly 1 (a hyperplane section).  So the cost of an order of the
56 x's is the number of blocked meetings.  I search orders (randomised greedy that prefers unblocked x's, then a
simulated-annealing polish) and report, per support, the best cost found; also the count-only floor (runs - 1, runs
from the claim's exact maxima: 4 natural, 8 best twist) and the final frame dimension of the gather order (what the
U-drop would then cost, instead of the 8 the claim charges separately).
Cases: the natural mirror (A = I) and the claim's exhaustive best twist (A = I + M, M from twistall.json)."""
import gzip, json, os, random, sys, time
T0 = time.time()
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(gzip.open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))))
U = list(d['ports']); v = d['v']; J = (1 << 9) - 1
supp = json.load(open(os.path.join(here, 'v1.json')))['supp']
def dot(a, b): return bin(a & b).count('1') & 1
def apply(A, x):
    y = 0
    for i in range(9):
        if x >> i & 1: y ^= A[i]
    return y
tw = json.load(open(os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../second-line-mirrored/all-twists/twistall.json'))))
M = tw['best']['M']
cases = {'natural': [1 << i for i in range(9)], 'best_twist': [(1 << i) ^ M[i] for i in range(9)]}
def reduce_basis(B):
    out = []
    for x in B:
        for y in out: x = min(x, x ^ y)
        if x: out.append(x)
    return out
def sect(F, g):
    """F & g^perp (F a reduced list of window vectors)"""
    bad = [f for f in F if dot(f, g)]
    if not bad: return F, 0
    p = bad[0]
    return reduce_basis([f ^ p if dot(f, g) else f for f in F if f != p]), 1
def run(order, G):
    F = []; cost = 0
    for b in order:
        F, c = sect(F, G[b]); cost += c
        F = reduce_basis(F + [U[b]])
    return cost, len(F)
def greedy(S, G, rnd):
    rest = list(S); F = []; cost = 0; order = []
    while rest:
        free = [b for b in rest if all(dot(f, G[b]) == 0 for f in F)]
        if free:
            # prefer the x that keeps the most remaining x's unblocked afterwards
            rnd.shuffle(free)
            def score(b):
                F2 = reduce_basis(F + [U[b]])
                return sum(1 for c in rest if c != b and all(dot(f, G[c]) == 0 for f in F2))
            b = max(free[:12], key=score)
        else:
            cands = rest[:]; rnd.shuffle(cands)
            def score2(b):
                F2, _ = sect(F, G[b]); F2 = reduce_basis(F2 + [U[b]])
                return sum(1 for c in rest if c != b and all(dot(f, G[c]) == 0 for f in F2))
            b = max(cands[:16], key=score2)
            F, c = sect(F, G[b]); cost += c
        F = reduce_basis(F + [U[b]]); order.append(b); rest.remove(b)
    return cost, order
def anneal(order, G, rnd, iters):
    cur = list(order); cc, _ = run(cur, G); best = (cc, list(cur))
    for it in range(iters):
        i, j = rnd.randrange(len(cur)), rnd.randrange(len(cur))
        if i == j: continue
        nxt = list(cur)
        if rnd.random() < 0.5: nxt[i], nxt[j] = nxt[j], nxt[i]
        else: x = nxt.pop(i); nxt.insert(j, x)
        nc, _ = run(nxt, G)
        temp = 0.6 * (1 - it / iters) + 0.05
        if nc <= cc or rnd.random() < pow(2.718, (cc - nc) / temp):
            cur, cc = nxt, nc
            if cc < best[0]: best = (cc, list(cur))
    return best
res = {}
budget = float(sys.argv[1]) if len(sys.argv) > 1 else 90.0
for name, A in cases.items():
    G = [J ^ apply(A, u) for u in U]          # g_b = j + A u_b  (window part of v_b)
    assert all(dot(U[b], G[b]) == 0 for b in range(v))
    rnd = random.Random(1)
    per = []
    for k, S in enumerate(supp):
        best = None
        t1 = time.time()
        while time.time() - t1 < budget / 16:
            c, o = greedy(S, G, rnd)
            if best is None or c < best[0]: best = (c, o)
        c2, o2 = anneal(best[1], G, rnd, 1500)
        if c2 < best[0]: best = (c2, o2)
        cost, fdim = run(best[1], G)
        assert cost == best[0]
        per.append(dict(cost=cost, final_dim=fdim, order=best[1]))
    res[name] = per
    print('%-11s lazy one-level drops per total (best order found): %s | sum %d per call per phase | final frame dims %s'
          % (name, [p['cost'] for p in per], sum(p['cost'] for p in per), [p['final_dim'] for p in per]))
json.dump(res, open(os.path.join(here, 'v3order.json'), 'w'))
print('secs %.1f' % (time.time() - T0))
