"""v-R1-recheck: independent encoding of the 2-layer relay relaxation (no closure), no gauge.
B = X Y^T over F_2, rank <= r.  layer bit z_S in {0,1}.  For S != T unknown: B[S,T] = 1 -> (z_T = 0 and z_S = 1).
B[S,S] = 1.  Also an optional symmetry-free sanity: rank r+1 static must be SAT.
usage: venv/bin/python -I indep_l2.py INSTANCE r layers(1|2) [timeout]"""
import sys, time, threading
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
from patterns import instance
from pysat.solvers import Glucose4
name, r, L = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
tmo = float(sys.argv[4]) if len(sys.argv) > 4 else 100
knows, _ = instance(name); v = len(knows)
n = [0]
def new():
    n[0] += 1; return n[0]
X = [[new() for _ in range(r)] for _ in range(v)]
Y = [[new() for _ in range(r)] for _ in range(v)]
ge = [[new() for _ in range(max(L - 1, 0))] for _ in range(v)]
cl = []
for S in range(v):
    for l in range(L - 2): cl.append([-ge[S][l + 1], ge[S][l]])
from patterns import max_clique_nonknow
K = max_clique_nonknow(knows, size=r)
assert all(not knows[a][b] for a in K for b in K if a != b)
for i, kk in enumerate(K):            # X_K = [I 0]: WLOG since B[K,K] is unitriangular in layer order
    for k in range(r): cl.append([X[kk][k]] if k == i else [-X[kk][k]])
def bentry(S, T):
    prods = []
    for k in range(r):
        p = new(); a, b = X[S][k], Y[T][k]
        cl.extend([[-p, a], [-p, b], [p, -a, -b]]); prods.append(p)
    acc = prods[0]
    for p in prods[1:]:
        q = new(); cl.extend([[-q, acc, p], [-q, -acc, -p], [q, -acc, p], [q, acc, -p]]); acc = q
    return acc
t0 = time.time()
for S in range(v):
    cl.append([bentry(S, S)])
    for T in range(v):
        if T != S and not knows[S][T]:
            b = bentry(S, T)
            if L == 1: cl.append([-b])
            else:   # layer(T) < layer(S): some threshold l with layer(S) >= l > layer(T)
                ws = []
                for l in range(L - 1):
                    w = new(); cl.extend([[-w, ge[S][l]], [-w, -ge[T][l]]]); ws.append(w)
                cl.append([-b] + ws)
s = Glucose4(bootstrap_with=cl)
tm = threading.Timer(tmo, s.interrupt); tm.start()
res = s.solve_limited(expect_interrupt=True); tm.cancel()
print('%s r=%d L=%d X_K=[I 0] |K|=%d, independent encoding (Glucose4): %s  vars %d clauses %d  %.1fs' % (name, r, L, len(K), {True: 'SAT', False: 'UNSAT', None: 'UNDECIDED'}[res], n[0], len(cl), time.time() - t0))
