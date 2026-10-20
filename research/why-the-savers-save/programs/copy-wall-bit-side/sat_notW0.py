"""sat_notW0.py (d-bitspan) -- does the bit pattern on P111 have an F_2 decoder with K totals whose read span W does
NOT contain the label space W0 = span{r(q)}?  No flat constraints at all (any read sets, any cost).
Model: X (K x v), Y (K x v) over F_2, B = Y^T X fits (B[T,T] = 1, B[S,T] = 0 on non-orthogonal live pairs);
z in F_2^v with z . x_j = 0 for all j (z in W^perp) and z . r(e_k) = 1 for some k (z not in W0^perp).
Gauge: X restricted to an 8-clique C is [I_8; 0] (GL_K acts on (X, Y) by (A X, A^-T Y); rank X|_C = 8).
UNSAT at K means: every decoder with <= K totals has W containing W0 ('case A').
usage: venv/bin/python -I sat_notW0.py <s2|word> K"""
import sys, time; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from e8geo import *
unit = sys.argv[1]; K = int(sys.argv[2])
P = {'s2': PUNCT_S2, 'word': PUNCT_WORD}[unit]
live = [i for i in range(NL) if i not in P]; v = len(live)
U = f2_labels()
C = cliques(live, 8)[0]; cpos = {live.index(T): i for i, T in enumerate(C)}
pool = IDPool(); cl = []
TRUE = pool.id('TRUE'); cl.append([TRUE])
X = [[None] * v for _ in range(K)]
for j in range(K):
    for t in range(v):
        if t in cpos: X[j][t] = TRUE if cpos[t] == j else -TRUE
        else: X[j][t] = pool.id(('x', j, t))
Y = [[(TRUE if cpos[s] == j else -TRUE) if s in cpos else pool.id(('y', j, s)) for s in range(v)] for j in range(K)]   # gauge Y|_C = [I_8; 0] too (as the hub's solve_asym)
def AND(a, b):
    if a == TRUE: return b
    if b == TRUE: return a
    if a == -TRUE or b == -TRUE: return -TRUE
    z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def XORc(lits, rhs):
    lits = [l for l in lits if l != -TRUE]
    k = sum(1 for l in lits if l == TRUE); lits = [l for l in lits if l != TRUE]
    rhs ^= k & 1
    if not lits:
        if rhs: cl.append([-TRUE])
        return
    acc = lits[0]
    for l in lits[1:]:
        z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
    cl.append([acc] if rhs else [-acc])
def XORv(lits):
    acc = lits[0]
    for l in lits[1:]:
        z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
    return acc
npairs = 0
for s in range(v):
    for t in range(v):
        if s == t or G[live[s], live[t]] != 0:
            XORc([AND(Y[j][s], X[j][t]) for j in range(K)], 1 if s == t else 0); npairs += 1
z = [pool.id(('z', t)) for t in range(v)]
for j in range(K):
    XORc([AND(z[t], X[j][t]) for t in range(v)], 0)
a = [XORv([z[t] for t in range(v) if U[live[t]][k]]) for k in range(9)]
cl.append(a)
print('unit %s v %d K %d gauge clique %s vars %d clauses %d fit pairs %d' % (unit, v, K, C, pool.top, len(cl), npairs), flush=True)
s = Cadical153(bootstrap_with=cl); t0 = time.time(); r = s.solve()
print('RESULT', 'SAT (a decoder whose read span misses W0)' if r else 'UNSAT (every decoder with K totals has W >= W0)', '%.1f s' % (time.time() - t0), flush=True)
if r:
    m = set(l for l in s.get_model() if l > 0)
    import json
    Xs = [[1 if (X[j][t] == TRUE or (X[j][t] not in (TRUE, -TRUE) and X[j][t] in m)) else 0 for t in range(v)] for j in range(K)]
    Ys = [[1 if (Y[j][t] == TRUE or (Y[j][t] not in (TRUE, -TRUE) and Y[j][t] in m)) else 0 for t in range(v)] for j in range(K)]
    json.dump(dict(unit=unit, live=live, X=Xs, Y=Ys), open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'out/notW0_%s_K%d.json') % (unit, K), 'w'))
    sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-floor-checker'))
    import copyfloor as cf
    Ub = cf.Unit('e8-bit-p111', 'Q', [list(map(int, L[i])) for i in live], p=2)
    terms = list(zip(Xs, Ys)); bad = Ub.fits(Ub.delivered(terms)); st, pr, rows = Ub.cost(terms)
    print('copyfloor check: fits %s, star %d, pairing %d, per-total star %s' % (not bad, st, pr, [rr['star'] for rr in rows]), flush=True)
