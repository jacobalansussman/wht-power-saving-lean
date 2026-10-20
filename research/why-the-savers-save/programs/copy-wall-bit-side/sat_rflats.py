"""sat_rflats.py (d-bitspan) -- off-pattern F_2 decoders of the BIT pattern on a punctured E8 unit, with every total's
reads confined to a root-orthogonal flat.
Total j has a class c_j in {7, 6, 5, 8}:
  7: reads inside l_j^perp            (E7 hyperplane, star cost <= 7), l_j any of the 120 lines;
  6: reads inside l_j^perp cap m_j^perp, l_j != m_j   (E6 or D6 6-flat, star cost <= 6);
  5: reads inside (l_j, m_j, n_j)^perp, three lines spanning rank 3 (cost <= 5);
  8: reads anything (cost 8).
Fit over F_2 (any rank K = number of totals): B = sum_j y_j x_j^T, B[T,T] = 1, B[S,T] = 0 when S != T are
non-orthogonal (both live).  B is otherwise free (orthogonal pairs) -- the hub's 'off-pattern' decoders included.
Symmetry breaking: totals of one class sorted by l_j; l_j < m_j (< n_j).
usage: venv/bin/python -I sat_rflats.py <unit: s2|word|120> <profile e.g. 6x8,7x2> [--noperm]"""
import sys, time, itertools; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from e8geo import *

def build(live, classes, punct_only=None):
    pool = IDPool(); cl = []
    v = len(live); K = len(classes)
    X = [[pool.id(('x', j, t)) for t in range(v)] for j in range(K)]
    Y = [[pool.id(('y', j, s)) for s in range(v)] for j in range(K)]
    sel = []
    cand = list(range(NL)) if punct_only is None else punct_only
    for j, c in enumerate(classes):
        nsel = {8: 0, 7: 1, 6: 2, 5: 3}[c]
        S = [[pool.id(('sel', j, a, l)) for l in cand] for a in range(nsel)]
        sel.append(S)
        for a in range(nsel):
            cl.append(S[a])                                         # at least one
            for p in range(len(cand)):
                for q in range(p + 1, len(cand)):
                    cl.append([-S[a][p], -S[a][q]])                 # at most one
            for t in range(v):                                      # reads only lines orthogonal to the chosen line
                T = live[t]
                cl.append([-X[j][t]] + [S[a][k] for k, l in enumerate(cand) if G[T, l] == 0 and l != T])
        # distinct, increasing selections l < m < n; and rank: for class 6 need l != m (any two distinct lines
        # are independent); for class 5 require rank 3 (forbid a third line inside span(l, m), i.e. the A2 third line)
        for a in range(1, nsel):
            for p in range(len(cand)):
                for q in range(p + 1):
                    cl.append([-S[a][p], -S[a - 1][q]] if q <= p else [])
        if nsel == 3:
            for p, q, r in itertools.combinations(range(len(cand)), 3):
                l, m, n = cand[p], cand[q], cand[r]
                M = np.array([L[l], L[m], L[n]])
                if np.linalg.matrix_rank(M.astype(float)) < 3:
                    cl.append([-S[0][p], -S[1][q], -S[2][r]])
    cl = [c for c in cl if c]
    # order totals of the same class by first selector
    for j in range(K - 1):
        if classes[j] == classes[j + 1] and classes[j] != 8:
            A, B_ = sel[j][0], sel[j + 1][0]
            for p in range(len(cand)):
                for q in range(p):                                  # forbid sel_{j+1} < sel_j
                    cl.append([-A[p], -B_[q]])
    # fit
    def AND(a, b):
        z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
    def XORc(lits, rhs):
        acc = lits[0]
        for l in lits[1:]:
            z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
        cl.append([acc] if rhs else [-acc])
    npairs = 0
    for s in range(v):
        for t in range(v):
            S_, T_ = live[s], live[t]
            if s == t or G[S_, T_] != 0:
                XORc([AND(Y[j][s], X[j][t]) for j in range(K)], 1 if s == t else 0); npairs += 1
    return pool, cl, X, Y, sel, cand, npairs

if __name__ == '__main__':
    unit = sys.argv[1]
    P = {'s2': PUNCT_S2, 'word': PUNCT_WORD, '120': []}[unit]
    live = [i for i in range(NL) if i not in P]
    classes = []
    for part in sys.argv[2].split(','):
        c, n = part.split('x'); classes += [int(c)] * int(n)
    classes.sort(reverse=True)
    t0 = time.time()
    pool, cl, X, Y, sel, cand, npairs = build(live, classes)
    print('unit %s v %d classes %s (star cost <= %d), vars %d clauses %d, fit pairs %d, built %.1f s'
          % (unit, len(live), classes, sum(classes), pool.top, len(cl), npairs, time.time() - t0), flush=True)
    s = Cadical153(bootstrap_with=cl); t1 = time.time(); r = s.solve()
    print('RESULT', 'SAT' if r else 'UNSAT', '%.1f s' % (time.time() - t1), flush=True)
    if r:
        m = set(l for l in s.get_model() if l > 0)
        import json
        terms = []
        for j in range(len(classes)):
            chosen = [[cand[k] for k in range(len(cand)) if S[k] in m] for S in sel[j]]
            R = [live[t] for t in range(len(live)) if X[j][t] in m]
            D = [live[t] for t in range(len(live)) if Y[j][t] in m]
            terms.append(dict(cls=classes[j], lines=chosen, read=R, served=D))
            print(' total', j, 'class', classes[j], 'flat lines', chosen, 'reads', len(R), 'serves', len(D), flush=True)
        json.dump(dict(unit=unit, live=live, terms=terms), open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/sat_%s_%s.json') % (unit, sys.argv[2].replace(',', '_')), 'w'))
