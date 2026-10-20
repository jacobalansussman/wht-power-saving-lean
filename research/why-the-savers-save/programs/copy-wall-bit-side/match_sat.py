"""match_sat.py (d-bitspan) -- the Hall/matching relaxation of the bit copy floor, field-free.
Every fitting decoder (ANY payload field, ANY rank) has B[C,C] = I_8 on every bit clique C (8 pairwise
non-orthogonal live lines), so the read vectors restricted to C span K^C; the restriction of x_j is supported on
C cap F_j, hence (Rado) the bipartite graph totals x lines(C), j ~ T iff T in F_j, has a matching saturating C.
SAT over the choice of flats: is there a family of K totals with classes c_j (7: inside l^perp; 6: inside
(l,m)^perp, l != m; 5: inside (l,m,n)^perp rank 3; 8: anything) such that EVERY clique has a saturating matching?
UNSAT = no decoder with that cost profile in the root-orthogonal flat family, over any field and at any rank.
usage: venv/bin/python -I match_sat.py <s2|word|120> <profile e.g. 6x8,7x2> [--cand punct]"""
import sys, time, itertools, json; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from e8geo import *

def build(live, classes, cand=None):
    pool = IDPool(); cl = []
    cand = list(range(NL)) if cand is None else cand
    K = len(classes)
    cov = {}
    sel = []
    for j, c in enumerate(classes):
        nsel = {8: 0, 7: 1, 6: 2, 5: 3}[c]
        S = [[pool.id(('sel', j, a, l)) for l in cand] for a in range(nsel)]
        sel.append(S)
        for a in range(nsel):
            enc = CardEnc.equals(lits=S[a], bound=1, vpool=pool, encoding=EncType.seqcounter); cl.extend(enc.clauses)
        for a in range(1, nsel):                                 # l < m < n
            for p in range(len(cand)):
                for q in range(p, len(cand)):
                    cl.append([-S[a][p], -S[a - 1][q]])
        if nsel == 3:
            for p, q, r in itertools.combinations(range(len(cand)), 3):
                M = np.array([L[cand[p]], L[cand[q]], L[cand[r]]], dtype=float)
                if np.linalg.matrix_rank(M) < 3: cl.append([-S[0][p], -S[1][q], -S[2][r]])
        for T in live:                                           # cov(j,T) -> T orthogonal to every chosen line
            z = pool.id(('cov', j, T)); cov[j, T] = z
            for a in range(nsel):
                cl.append([-z] + [S[a][k] for k, l in enumerate(cand) if G[T, l] == 0 and l != T])
    for j in range(K - 1):                                       # same class: sorted by first selector
        if classes[j] == classes[j + 1] and classes[j] != 8:
            A, B_ = sel[j][0], sel[j + 1][0]
            for p in range(len(cand)):
                for q in range(p):
                    cl.append([-A[p], -B_[q]])
    C = cliques(live, 8)
    for ci, c in enumerate(C):
        m = {}
        for T in c:
            for j in range(K):
                m[T, j] = pool.id(('m', ci, T, j)); cl.append([-m[T, j], cov[j, T]])
            cl.append([m[T, j] for j in range(K)])
        for j in range(K):
            lits = [m[T, j] for T in c]
            for a in range(8):
                for b in range(a + 1, 8): cl.append([-lits[a], -lits[b]])
    # every live line covered (diagonal) -- implied by cliques if every line is in a clique; add anyway
    for T in live:
        cl.append([cov[j, T] for j in range(K)])
    return pool, cl, sel, cand, len(C)

if __name__ == '__main__':
    unit = sys.argv[1]
    P = {'s2': PUNCT_S2, 'word': PUNCT_WORD, '120': []}[unit]
    live = [i for i in range(NL) if i not in P]
    classes = []
    for part in sys.argv[2].split(','):
        c, n = part.split('x'); classes += [int(c)] * int(n)
    classes.sort(reverse=True)
    cand = P if '--cand' in sys.argv and sys.argv[sys.argv.index('--cand') + 1] == 'punct' else None
    t0 = time.time()
    pool, cl, sel, cand, ncl = build(live, classes, cand)
    print('unit %s v %d cliques %d classes %s (cost <= %d) vars %d clauses %d built %.1f s'
          % (unit, len(live), ncl, classes, sum(classes), pool.top, len(cl), time.time() - t0), flush=True)
    s = Cadical153(bootstrap_with=cl); t1 = time.time(); r = s.solve()
    print('RESULT', 'SAT' if r else 'UNSAT', '%.1f s' % (time.time() - t1), flush=True)
    if r:
        mset = set(l for l in s.get_model() if l > 0)
        fam = []
        for j in range(len(classes)):
            ch = [[cand[k] for k in range(len(cand)) if S[k] in mset][0] for S in sel[j]]
            fam.append(dict(cls=classes[j], lines=ch))
            print(' total', j, 'class', classes[j], 'flat = perp of lines', ch, 'types', [(a, b, int(G[a, b])) for a, b in itertools.combinations(ch, 2)], flush=True)
        json.dump(dict(unit=unit, fam=fam), open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/match_%s_%s.json') % (unit, sys.argv[2].replace(',', '_')), 'w'))
