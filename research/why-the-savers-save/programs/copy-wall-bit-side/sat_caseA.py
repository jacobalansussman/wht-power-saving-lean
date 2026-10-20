"""sat_caseA.py (d-bitspan) -- CASE A of the off-pattern bit decoders on P111: decoders whose delivered matrix is still
B0 but whose read vectors leave the label space W0 = span of the 9 label functionals r(q), q in F_2^9.
Every decoder with K = 10 totals and B = B0 has span{x_j} = W0 + <w> for one extra vector w in F_2^v; so
x_j = r(q_j) + c_j w.  Validity: the rows (q_j, c_j) span F_2^9 x {0} (then r(e_k) = sum of x_j, so W0 is inside
the span and served vectors exist).  The hub's label-decoder SAT (dec9sat/run62, our graphsat) covers only c_j = 0.
Read confinement as in sat_rflats.py (class 7: l^perp; 6: (l,m)^perp; 5: (l,m,n)^perp; 8: free).
usage: venv/bin/python -I sat_caseA.py <s2|word> <profile e.g. 6x8,7x2> [--hyp punct]   (writes out/caseA_*.json if SAT)"""
import sys, time, itertools, json; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from e8geo import *

def build(live, classes, U):
    pool = IDPool(); cl = []
    v = len(live); K = len(classes)
    q = [[pool.id(('q', j, k)) for k in range(9)] for j in range(K)]
    c = [pool.id(('c', j)) for j in range(K)]
    w = [pool.id(('w', t)) for t in range(v)]
    def AND(a, b):
        z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
    def XOR(lits):
        if not lits: return None
        acc = lits[0]
        for l in lits[1:]:
            z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
        return acc
    X = []
    for j in range(K):
        row = []
        for t in range(v):
            T = live[t]
            row.append(XOR([q[j][k] for k in range(9) if U[T][k]] + [AND(c[j], w[t])]))
        X.append(row)
    sel = []
    for j, cc in enumerate(classes):
        nsel = {8: 0, 7: 1, 6: 2, 5: 3}[cc]
        S = [[pool.id(('sel', j, a, l)) for l in range(NL)] for a in range(nsel)]
        sel.append(S)
        for a in range(nsel):
            cl.extend(CardEnc.equals(lits=S[a], bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
            for t in range(v):
                T = live[t]
                cl.append([-X[j][t]] + [S[a][l] for l in range(NL) if G[T, l] == 0 and l != T])
        for a in range(1, nsel):
            for p in range(NL):
                for r in range(p, NL): cl.append([-S[a][p], -S[a - 1][r]])
        if nsel == 3:
            for p, r, s in itertools.combinations(range(NL), 3):
                if np.linalg.matrix_rank(np.array([L[p], L[r], L[s]], float)) < 3: cl.append([-S[0][p], -S[1][r], -S[2][s]])
    # validity: N (9 x K) with N A = [I_9 | 0], A rows (q_j, c_j)
    Nv = [[pool.id(('N', k, j)) for j in range(K)] for k in range(9)]
    for k in range(9):
        for i in range(10):
            col = [q[j][i] for j in range(K)] if i < 9 else c
            s = XOR([AND(Nv[k][j], col[j]) for j in range(K)])
            cl.append([s] if i == k else [-s])
    # symmetry: same-class totals sorted by first selector
    for j in range(K - 1):
        if classes[j] == classes[j + 1] and classes[j] != 8:
            A_, B_ = sel[j][0], sel[j + 1][0]
            for p in range(NL):
                for r in range(p): cl.append([-A_[p], -B_[r]])
    return pool, cl, q, c, w, X, sel

def decode(live, U, qs, cs, wv):
    """read vectors x_j and served vectors y_j with sum_j y_j x_j^T = B0|live (solve over F_2)."""
    v = len(live)
    X = [[(sum(U[live[t]][k] * qs[j][k] for k in range(9)) + cs[j] * wv[t]) % 2 for t in range(v)] for j in range(len(qs))]
    B0l = [[int((1 + G[live[s], live[t]]) % 2) for t in range(v)] for s in range(v)]
    K = len(X)
    # Gaussian elimination: express each row of B0l in the basis X (X may be dependent: use a spanning solve)
    rows = [(int(''.join(map(str, X[j][::-1])), 2), 1 << j) for j in range(K)]
    basis = {}
    for val, comb in rows:
        while val:
            h = val.bit_length() - 1
            if h in basis: bv, bc = basis[h]; val ^= bv; comb ^= bc
            else: basis[h] = (val, comb); break
    Y = [[0] * v for _ in range(K)]
    for s in range(v):
        val = int(''.join(map(str, B0l[s][::-1])), 2); comb = 0
        while val:
            h = val.bit_length() - 1
            if h not in basis: return X, None
            bv, bc = basis[h]; val ^= bv; comb ^= bc
        for j in range(K):
            if comb >> j & 1: Y[j][s] = 1
    return X, Y

if __name__ == '__main__':
    unit = sys.argv[1]
    P = {'s2': PUNCT_S2, 'word': PUNCT_WORD}[unit]
    live = [i for i in range(NL) if i not in P]
    classes = []
    for part in sys.argv[2].split(','):
        cc, n = part.split('x'); classes += [int(cc)] * int(n)
    classes.sort(reverse=True)
    U = f2_labels()
    t0 = time.time()
    pool, cl, q, c, w, X, sel = build(live, classes, U)
    print('CASE A unit %s v %d classes %s (cost <= %d) vars %d clauses %d built %.1f s' % (unit, len(live), classes, sum(classes), pool.top, len(cl), time.time() - t0), flush=True)
    s = Cadical153(bootstrap_with=cl); t1 = time.time(); r = s.solve()
    print('RESULT', 'SAT' if r else 'UNSAT', '%.1f s' % (time.time() - t1), flush=True)
    if r:
        m = set(l for l in s.get_model() if l > 0)
        qs = [[int(q[j][k] in m) for k in range(9)] for j in range(len(classes))]
        cs = [int(c[j] in m) for j in range(len(classes))]
        wv = [int(w[t] in m) for t in range(len(live))]
        Xr, Y = decode(live, U, qs, cs, wv)
        sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker'))
        import copyfloor as cf
        Ub = cf.Unit('e8-bit-p111', 'Q', [list(map(int, L[i])) for i in live], p=2)
        terms = list(zip(Xr, Y))
        bad = Ub.fits(Ub.delivered(terms))
        st, pr, rows = Ub.cost(terms)
        print('copyfloor: fits %s, star %d, pairing %d, per-total star %s, c_j %s, |w| %d' % (not bad, st, pr, [rr['star'] for rr in rows], cs, sum(wv)), flush=True)
        json.dump(dict(unit=unit, live=live, q=qs, c=cs, w=wv, X=Xr, Y=Y, star=st, pair=pr),
                  open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/caseA_%s_%s.json') % (unit, sys.argv[2].replace(',', '_')), 'w'))
