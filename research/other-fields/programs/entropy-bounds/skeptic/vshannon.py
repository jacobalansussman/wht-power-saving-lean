# vshannon.py -- v-nc-gain's own Shannon (polymatroid) LP for min-cost network coding (independent of nclp.py).
# Ground set: sources Y_a (rate 1, independent) and arcs X_e. Rows: elemental Shannon inequalities; source
# independence; per-node closure h(In(v)) = h(In(v) + Out(v) + Dec(v)) (In = sources at v + arcs into v);
# objective min sum_e h(X_e). Optional Ingleton rows (all 4-tuples of single elements with empty context,
# A,B ; C,D) as a check that they do not change calibration values.
import itertools, sys
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

def shannon_lp(nodes, arcs, sessions, ingleton=False):
    """sessions: list of (src node, [sink nodes]). returns LP value."""
    ns = len(sessions); na = len(arcs); n = ns + na; N = 1 << n
    Y = lambda a: 1 << a; X = lambda e: 1 << (ns + e)
    rows = []   # list of dict mask->coef meaning sum coef*h(mask) >= 0
    full = N - 1
    for i in range(n): rows.append({full: 1, full ^ (1 << i): -1})
    for i, j in itertools.combinations(range(n), 2):
        rest = [k for k in range(n) if k not in (i, j)]
        for r in range(len(rest) + 1):
            for K in itertools.combinations(rest, r):
                k = sum(1 << x for x in K)
                d = {}
                for msk, cf in ((k | 1 << i, 1), (k | 1 << j, 1), (k | 1 << i | 1 << j, -1), (k, -1)):
                    if msk: d[msk] = d.get(msk, 0) + cf
                rows.append(d)
    if ingleton:
        for A, B, C, D in itertools.permutations(range(n), 4):
            if A > B or C > D: continue
            a, b, c, dd = 1 << A, 1 << B, 1 << C, 1 << D
            # Ingleton: h(A)+h(B)+h(CD)+h(ABC)+h(ABD) <= h(AB)+h(AC)+h(AD)+h(BC)+h(BD)
            row = {}
            for msk, cf in ((a, -1), (b, -1), (c | dd, -1), (a | b | c, -1), (a | b | dd, -1),
                            (a | b, 1), (a | c, 1), (a | dd, 1), (b | c, 1), (b | dd, 1)):
                row[msk] = row.get(msk, 0) + cf
            rows.append(row)
    eqs = []   # (dict, rhs)
    allY = sum(Y(a) for a in range(ns))
    for a in range(ns): eqs.append(({Y(a): 1}, 1.0))
    eqs.append(({allY: 1}, float(ns)))
    for v in nodes:
        In = sum(Y(a) for a, (s, _) in enumerate(sessions) if s == v) | sum(X(e) for e, (p, q) in enumerate(arcs) if q == v)
        Out = sum(X(e) for e, (p, q) in enumerate(arcs) if p == v)
        Dec = sum(Y(a) for a, (_, ts) in enumerate(sessions) if v in ts)
        Bg = In | Out | Dec
        if Bg == In: continue
        if In == 0:
            eqs.append(({Bg: 1}, 0.0))
        else:
            eqs.append(({Bg: 1, In: -1}, 0.0))
    def mat(lst):
        I, J, V = [], [], []
        for r, d in enumerate(lst):
            for msk, cf in d.items():
                if cf: I.append(r); J.append(msk - 1); V.append(cf)
        return sp.csr_matrix((V, (I, J)), shape=(len(lst), N - 1))
    Aub = -mat(rows); bub = np.zeros(len(rows))
    Aeq = mat([d for d, _ in eqs]); beq = np.array([r for _, r in eqs])
    c = np.zeros(N - 1)
    for e in range(na): c[X(e) - 1] = 1.0
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
    return res.fun

if __name__ == '__main__':
    print("path 0-1-2-3 (directed), one session: expect 3 =", shannon_lp([0, 1, 2, 3], [(0, 1), (1, 2), (2, 3)], [(0, [3])]))
    bf = [('s1', 't1'), ('s1', 'm'), ('s2', 'm'), ('s2', 't2'), ('m', 'n'), ('n', 't1'), ('n', 't2')]
    nodes = ['s1', 's2', 'm', 'n', 't1', 't2']
    print("directed butterfly, 2 sources multicast to 2 sinks: expect 7 (routing 8) =", shannon_lp(nodes, bf, [('s1', ['t1', 't2']), ('s2', ['t1', 't2'])]))
    print("  same with Ingleton rows: ", shannon_lp(nodes, bf, [('s1', ['t1', 't2']), ('s2', ['t1', 't2'])], ingleton=True))
    print("butterfly as two unicasts s1->t2, s2->t1 (side arcs s1->t1, s2->t2 useless): routing 6 =", shannon_lp(nodes, bf, [('s1', ['t2']), ('s2', ['t1'])]))
    comb = [('s', 'a'), ('s', 'b'), ('s', 'c'), ('a', 'tab'), ('b', 'tab'), ('a', 'tac'), ('c', 'tac'), ('b', 'tbc'), ('c', 'tbc')]
    cn = ['s', 'a', 'b', 'c', 'tab', 'tac', 'tbc']
    print("combination network C(3,2), 2 symbols multicast to 3 sinks: expect 9 (routing 10) =",
          shannon_lp(cn, comb, [('s', ['tab', 'tac', 'tbc']), ('s', ['tab', 'tac', 'tbc'])]))
    print("one edge both ways, sessions u->v, v->u: expect 2 =", shannon_lp(['u', 'v'], [('u', 'v'), ('v', 'u')], [('u', ['v']), ('v', ['u'])]))
