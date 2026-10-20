"""probe_puncture.py (g-break, H14/H5) -- linear-total star floor of E8 minus a set Z of lines.
Linear totals: c in Q^8, read set R(c) = {T live : <w_T, c> != 0}, star = F_2 span dim of u(R).  Any basis of
c's (with its dual basis as served vectors) gives a fitting decoder B = Gram on the live ports, so the linear
star floor = min-weight basis (matroid greedy) over candidate c's.  Candidates: all c in {-2..2}^8 up to sign
whose orthogonal complement contains >= 7 independent roots (otherwise R(c) is big) -- we enumerate c as the
normals of hyperplanes spanned by roots, i.e. c orthogonal to 7 independent roots, plus coordinate vectors."""
import itertools, sys, time
sys.dont_write_bytecode = True
from fractions import Fraction
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), ".")); from common import *

def rank_q(rows):
    return cf.rank(rows, 'Q')

def normals():
    # c = e_a, and c = root directions (lines themselves), and c = sums of two orthogonal roots etc.
    cs = set()
    for a in range(8):
        cs.add(tuple(1 if i == a else 0 for i in range(8)))
    for l in L:
        g = max(abs(x) for x in l)
        cs.add(tuple(x // (2 if g == 2 else 1) for x in l))
    # normals of root hyperplanes: pairs/triples of orthogonal roots summed
    for i in range(120):
        for j in range(i + 1, 120):
            if G[i][j] == 0:
                for s in (1, -1):
                    c = tuple(x + s * y for x, y in zip(L[i], L[j]))
                    cs.add(c)
    return [list(c) for c in cs]

def linear_floor(live, cands):
    w = []
    for c in cands:
        R = [T for T in live if sum(x * y for x, y in zip(L[T], c)) != 0]
        w.append((f2dim([UI[T] for T in R]), c, len(R)))
    w.sort(key=lambda t: t[0])
    chosen, rows = [], []
    for wt, c, n in w:
        if rank_q(rows + [c]) > len(rows):
            rows.append(c); chosen.append((wt, n))
        if len(rows) == 8:
            break
    return sum(x for x, _ in chosen), chosen, w

if __name__ == '__main__':
    t0 = time.time()
    cands = normals()
    print('candidate normals', len(cands))
    octads = []
    # find octads (8 mutually orthogonal lines) containing structure: coordinate octad = lines 0..7
    O1 = list(range(8))
    # second octads disjoint from O1 among sign lines: greedy clique search
    orth = [[G[i][j] == 0 and i != j for j in range(120)] for i in range(120)]
    def cliques(cands_, size, cur, out, limit):
        if len(out) >= limit: return
        if len(cur) == size:
            out.append(list(cur)); return
        for k, x in enumerate(cands_):
            if all(orth[x][y] for y in cur):
                cliques(cands_[k + 1:], size, cur + [x], out, limit)
    oct2 = []
    cliques(list(range(8, 120)), 8, [], oct2, 400)
    print('second octads found (first 400):', len(oct2))
    res = {}
    live = [i for i in range(120) if i not in O1]
    s, ch, _ = linear_floor(live, cands)
    print('P112 (E8 minus coordinate octad): linear star floor %d, weights %s' % (s, [c[0] for c in ch]))
    s, ch, _ = linear_floor(list(range(120)), cands)
    print('E8 (P120): linear star floor %d, weights %s' % (s, [c[0] for c in ch]))
    best = None
    seen = 0
    for O2 in oct2:
        live = [i for i in range(120) if i not in O1 and i not in O2]
        s, ch, _ = linear_floor(live, cands)
        seen += 1
        if best is None or s < best[0]:
            best = (s, O2, [c[0] for c in ch])
        if seen >= 60: break
    print('P104 (two octads removed), best of %d: linear star floor %d, weights %s, O2 %s' % (seen, best[0], best[2], best[1]))
    print('time %.1f s' % (time.time() - t0))
