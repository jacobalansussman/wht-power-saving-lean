"""sample_caseB.py (d-bitspan) -- sample F_2 decoders with K = 10 totals whose read span W misses the label space W0
(case B), and price each W exactly: the star cost of the best decoder with read span W is the minimum-weight basis
of W (weight of x = Q-rank of the lines in supp x; matroid greedy over all 1023 nonzero vectors of W).
Diversified by random unit assumptions on read variables.  usage: venv/bin/python -I sample_caseB.py <s2|word> N seed"""
import sys, random, time, collections; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import numpy as np
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from e8geo import *
unit, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
P = {'s2': PUNCT_S2, 'word': PUNCT_WORD}[unit]
live = [i for i in range(NL) if i not in P]; v = len(live); K = 10
U = f2_labels()
Lf = L[live].astype(float)
def wt(mask):
    idx = [t for t in range(v) if mask >> t & 1]
    return int(np.linalg.matrix_rank(Lf[idx])) if idx else 0
def minbasis(rows):
    vecs = []
    for c in range(1, 1 << len(rows)):
        x = 0
        for j in range(len(rows)):
            if c >> j & 1: x ^= rows[j]
        if x: vecs.append(x)
    vecs = list(set(vecs))
    W = sorted(((wt(x), x) for x in vecs))
    basis = {}; chosen = []
    for w_, x in W:
        y = x
        while y:
            h = y.bit_length() - 1
            if h in basis: y ^= basis[h]
            else: basis[h] = y; chosen.append(w_); break
    return sum(chosen), sorted(chosen), collections.Counter(w_ for w_, x in W)
# --- model (as sat_notW0.py, with gauge)
C = cliques(live, 8)[0]; cpos = {live.index(T): i for i, T in enumerate(C)}
pool = IDPool(); cl = []
TRUE = pool.id('TRUE'); cl.append([TRUE])
X = [[(TRUE if cpos[t] == j else -TRUE) if t in cpos else pool.id(('x', j, t)) for t in range(v)] for j in range(K)]
Y = [[(TRUE if cpos[s] == j else -TRUE) if s in cpos else pool.id(('y', j, s)) for s in range(v)] for j in range(K)]
def AND(a, b):
    if a == TRUE: return b
    if b == TRUE: return a
    if a == -TRUE or b == -TRUE: return -TRUE
    z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def XORc(lits, rhs):
    lits = [l for l in lits if l != -TRUE]; k = sum(1 for l in lits if l == TRUE); lits = [l for l in lits if l != TRUE]; rhs ^= k & 1
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
for s in range(v):
    for t in range(v):
        if s == t or G[live[s], live[t]] != 0:
            XORc([AND(Y[j][s], X[j][t]) for j in range(K)], 1 if s == t else 0)
z = [pool.id(('z', t)) for t in range(v)]
for j in range(K): XORc([AND(z[t], X[j][t]) for t in range(v)], 0)
cl.append([XORv([z[t] for t in range(v) if U[live[t]][k]]) for k in range(9)])
sv = Cadical153(bootstrap_with=cl)
rng = random.Random(seed)
free_x = [X[j][t] for j in range(K) for t in range(v) if t not in cpos]
best = None; hist = collections.Counter(); t0 = time.time()
for it in range(N):
    assum = [lit if rng.random() < 0.5 else -lit for lit in rng.sample(free_x, 6)]
    if not sv.solve(assumptions=assum): hist['unsat-assump'] += 1; continue
    m = set(l for l in sv.get_model() if l > 0)
    rows = []
    for j in range(K):
        x = 0
        for t in range(v):
            if X[j][t] == TRUE or (X[j][t] != -TRUE and X[j][t] in m): x |= 1 << t
        rows.append(x)
    cost, ws, dist = minbasis(rows)
    hist[cost] += 1
    if best is None or cost < best[0]: best = (cost, ws, sorted(dist.items()))
    print('sample %d: min-weight-basis star cost %d  weights %s  (%.0f s)' % (it, cost, ws, time.time() - t0), flush=True)
print('SUMMARY unit %s samples %d cost histogram %s best %s' % (unit, N, sorted(hist.items(), key=str), best), flush=True)
