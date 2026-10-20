"""kfree.py (d-relay) -- g8 diagnostic: how many freed (relayed) pairs does a rank-r F_2 decoder need?
B = X Y^T, B[S,S] = 1, B[S,T] free when S knows T, and B[S,T] != 0 on at most k NON-known cells (the freed pairs,
in any order; 2- and 3-cycles of freed pairs are forbidden, so for k <= 3 the freed pairs are acyclic = realisable by
an order, and B[K,K] is unitriangular: the gauge X_K = I is without loss of generality).  No closure: this is the
decoder-level count, i.e. the most a PAID relay (T descends to a frame inside a_S^perp, priced d(f,g)) could use.
With 'cols' as 4th argument the bound is on the number of distinct RELAYED PORTS T (columns with a freed cell):
each relayed T needs its own paid descent (rank >= 1, 2 per call per rank by hub C.7), whatever the number of readers.
For k >= 4 every directed cycle of freed cells inside K is forbidden explicitly (keeps the gauge valid); longer cycles
outside K are not forbidden (a relaxation: UNSAT stays a proof).
usage: venv/bin/python -I kfree.py INSTANCE r k [cols]"""
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from patterns import instance, max_clique_nonknow
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
name, r, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
knows, G = instance(name); v = len(knows)
K = max_clique_nonknow(knows, size=r)
pool = IDPool(); cl = []
X = [[pool.id(('x', S, j)) for j in range(r)] for S in range(v)]
Y = [[pool.id(('y', S, j)) for j in range(r)] for S in range(v)]
for i, kk in enumerate(K):
    for j in range(r): X[kk][j] = 'T' if j == i else 'F'
def AND(a, b):
    if a == 'F' or b == 'F': return 'F'
    if a == 'T': return b
    if b == 'T': return a
    z = pool.id(('and', a, b)); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def XOR(ls):
    ls = [l for l in ls if l != 'F']
    c = ls.count('T') % 2; ls = [l for l in ls if l != 'T']
    if not ls: return 'T' if c else 'F'
    acc = ls[0]
    for l in ls[1:]:
        z = pool.id(('xor', acc, l)); cl.extend([[-z, acc, l], [-z, -acc, -l], [z, -acc, l], [z, acc, -l]]); acc = z
    if c:
        z = pool.id(('not', acc)); cl.extend([[z, acc], [-z, -acc]]); acc = z
    return acc
def Bv(S, T): return XOR([AND(X[S][j], Y[T][j]) for j in range(r)])
t0 = time.time()
for S in range(v):
    d = Bv(S, S)
    if d == 'F': cl.append([])
    elif d != 'T': cl.append([d])
free = {}
for S in range(v):
    for T in range(v):
        if S != T and not knows[S][T]:
            b = Bv(S, T)
            if b == 'F': continue
            if b == 'T': b = pool.id(('const1',)); cl.append([b])
            free[S, T] = b
for (S, T), b in free.items():
    if (T, S) in free and S < T: cl.append([-b, -free[T, S]])
if k >= 3:
    for (S, T), b in free.items():
        for U in range(v):
            if (T, U) in free and (U, S) in free and S < T and S < U: cl.append([-b, -free[T, U], -free[U, S]])
import itertools
if k >= 4:
    for l in range(4, len(K) + 1):
        for sub in itertools.combinations(K, l):
            first = sub[0]
            for perm in itertools.permutations(sub[1:]):
                cyc = (first,) + perm
                lits = [free.get((cyc[i], cyc[(i + 1) % l])) for i in range(l)]
                if all(lits): cl.append([-x for x in lits])
cols = len(sys.argv) > 4 and sys.argv[4] == 'cols'
if cols:
    cv = {}
    for (S, T), b in free.items():
        if T not in cv: cv[T] = pool.id(('col', T))
        cl.append([-b, cv[T]])
    card = CardEnc.atmost(lits=list(cv.values()), bound=k, vpool=pool, encoding=EncType.seqcounter)
else:
    card = CardEnc.atmost(lits=list(free.values()), bound=k, vpool=pool, encoding=EncType.seqcounter)
cl.extend(card.clauses)
s = Cadical153(bootstrap_with=cl); res = s.solve()
extra = ''
if res:
    m = set(l for l in s.get_model() if l > 0)
    used = [(S, T) for (S, T), b in free.items() if b in m]
    extra = 'freed pairs %s' % used
fmt_ = ('%s r=%d at most %d RELAYED PORTS:' if cols else '%s r=%d at most %d freed non-known cells:') + ' %s  vars %d clauses %d  %.1fs %s'
print(fmt_ % (name, r, k, 'SAT' if res else 'UNSAT', pool.top, len(cl), time.time() - t0, extra), flush=True)
