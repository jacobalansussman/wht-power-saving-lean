"""k9_indep.py (v-g-9hyp-recheck) -- independent encoding of 'every F_2 decoder with K totals has W >= W0' on P111.
Differences from d-bitspan's sat_notW0.py: own E8 build; W0 from a basis of the F_2 row space of B0|live (not f2_labels);
gauge on the LAST 8-clique found (not the first); Tseitin via pysat's own XOR-free encoding of products; solver Glucose4.
usage: python -I k9_indep.py <s2|word> K"""
import sys, itertools, time
import numpy as np
from pysat.solvers import Glucose4
from pysat.formula import IDPool
P_ = {'s2': [39, 52, 57, 71, 76, 78, 103, 109, 114], 'word': [8, 17, 39, 51, 53, 58, 63, 72, 114]}[sys.argv[1]]; K = int(sys.argv[2])
gen = np.array([[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]])
code = {tuple((np.array(c) @ gen) % 2) for c in itertools.product([0,1], repeat=4)}
w4 = sorted(c for c in code if sum(c) == 4)
Ls = [tuple(2*(np.arange(8)==a).astype(int)) for a in range(8)]
for c in w4:
    idx = [i for i in range(8) if c[i]]
    for sg in itertools.product([1,-1], repeat=3):
        x = np.zeros(8, dtype=int); x[idx[0]] = 1
        for k in range(3): x[idx[k+1]] = sg[k]
        Ls.append(tuple(x))
L = np.array(Ls, dtype=np.int64); G = (L @ L.T) // 2
live = [i for i in range(120) if i not in P_]; v = len(live)
NO = [[s != t and G[live[s], live[t]] != 0 for t in range(v)] for s in range(v)]
# all 8-cliques (pairwise non-orthogonal), take the last
cl8 = []
def rec(cur, cand):
    if len(cur) == 8: cl8.append(tuple(cur)); return
    for t in cand:
        rec(cur + [t], [u for u in cand if u > t and NO[t][u]])
rec([], list(range(v)))
C = cl8[-1]; cpos = {t: i for i, t in enumerate(C)}
B0 = ((1 + G) % 2)[np.ix_(live, live)]
basis = []; brows = []
for r in B0:
    x = int(''.join(str(int(b)) for b in r[::-1]), 2); y = x
    for b in basis: y = min(y, y ^ b)
    if y: basis.append(y); basis.sort(reverse=True); brows.append(r)
assert len(brows) == 9
pool = IDPool(); cls = []
T_ = pool.id('T'); cls.append([T_])
def lit_const(b): return T_ if b else -T_
X = [[lit_const(cpos[t] == j) if t in cpos else pool.id(('x', j, t)) for t in range(v)] for j in range(K)]
Y = [[lit_const(cpos[s] == j) if s in cpos else pool.id(('y', j, s)) for s in range(v)] for j in range(K)]
def AND(a, b):
    if a == -T_ or b == -T_: return -T_
    if a == T_: return b
    if b == T_: return a
    z = pool.id(('a', a, b)); cls.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def parity(lits):
    """literal equal to the XOR of lits (constants folded)."""
    par = 0; rest = []
    for l in lits:
        if l == T_: par ^= 1
        elif l == -T_: pass
        else: rest.append(l)
    if not rest: return T_ if par else -T_
    # balanced tree
    while len(rest) > 1:
        nxt = []
        for i in range(0, len(rest) - 1, 2):
            a, b = rest[i], rest[i + 1]; z = pool.id(('p', a, b))
            cls.extend([[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]); nxt.append(z)
        if len(rest) % 2: nxt.append(rest[-1])
        rest = nxt
    return -rest[0] if par else rest[0]
for s in range(v):
    for t in range(v):
        if s == t or NO[s][t]:
            p = parity([AND(Y[j][s], X[j][t]) for j in range(K)])
            cls.append([p] if s == t else [-p])
z = [pool.id(('z', t)) for t in range(v)]
for j in range(K):
    cls.append([-parity([AND(z[t], X[j][t]) for t in range(v)])])
cls.append([parity([z[t] for t in range(v) if r[t]]) for r in brows])
print('k9_indep unit %s K %d gauge clique (last of %d) %s vars %d clauses %d' % (sys.argv[1], K, len(cl8), C, pool.top, len(cls)), flush=True)
s = Glucose4(bootstrap_with=cls); t0 = time.time(); r = s.solve()
print('RESULT', 'SAT (W misses W0)' if r else 'UNSAT (W >= W0 for every decoder with K totals)', '%.1f s' % (time.time() - t0), flush=True)
