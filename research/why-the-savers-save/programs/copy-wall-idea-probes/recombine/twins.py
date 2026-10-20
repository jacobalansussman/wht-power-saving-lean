"""twins.py (g-recombine) -- H1/H13: do two calls that share one frame space (twin ports, same labels) need less
than twice the copies?  Small cases inside E8 (bit side, Q frames): pick line sets whose orthogonality graph is a
5-cycle (the classic index-coding case where vector codes beat scalar ones: minrank 3, rate 5/2) and other small
patterns; compute the exact F_2 star floor for t = 1 and for t = 2 twins with copyfloor's branch and bound.
usage: python3 -I twins.py"""
import sys, itertools, time
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
import copyfloor as cf
L, G = cf.e8_lines()
def find_cycle(k):
    O = lambda a, b: G[a][b] == 0
    for a in [0]:
        for b in range(120):
            if b == a or not O(a, b): continue
            for c in range(120):
                if c in (a, b) or not O(b, c) or O(a, c): continue
                for d in range(120):
                    if d in (a, b, c) or not O(c, d) or O(a, d) or O(b, d): continue
                    for e in range(120):
                        if e in (a, b, c, d) or not O(d, e) or not O(e, a) or O(b, e) or O(c, e): continue
                        return [a, b, c, d, e]
for k in (5,):
    cyc = find_cycle(k)
    print('C%d in E8 lines:' % k, cyc, 'Gram:', [[G[a][b] for b in cyc] for a in cyc])
    labels = [L[i] for i in cyc]
    print('  Q-rank of the labels:', cf.rank(labels, 'Q'))
    for t in (1, 2):
        U = cf.Unit('C%d-x%d' % (k, t), 'Q', [a for a in labels for _ in range(t)], p=2)
        t0 = time.time()
        val, fam, nodes = cf.star_floor_f2(U)
        print('  t = %d: v %d, star floor %d (per copy %.2f), nodes %d, %.1fs' % (t, U.v, val, val / t, nodes, time.time() - t0))
        for r in fam:
            R = [T for T in range(U.v) if r >> T & 1]
            print('     read', R, 'star', U.star(R))
