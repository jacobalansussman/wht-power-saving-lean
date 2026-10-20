"""sat9.py (v-g-9hyp-model) -- independent re-encoding of 'every F_2 decoder with K totals on P111 has read span W
containing W0'.  Decoder: X (K x v reads), Y (K x v serves) over F_2, B = Y^T X with B[S,S] = 1 and B[S,T] = 0
for S != T non-orthogonal (live).  Seek z in W^perp \\ W0^perp.  Gauge on an 8-clique C chosen by seed (different
from d-bitspan's first clique): X|_C = Y|_C = [I_8; 0] (WLOG: A in GL_K with A X|_C = [I;0], then
A = [[I, A12],[0, I]] with A12 = Y_bot|_C^T clears Y_bot|_C, and B|_{CxC} = I forces Y_top|_C = I).
XORs encoded by balanced Tseitin trees.  usage: python -I sat9.py <s2|word> K seed [solver]"""
import sys, time, random; sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
import geo9
from pysat.solvers import Solver
name, K, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
solver = sys.argv[4] if len(sys.argv) > 4 else 'cadical153'
P, live, W0b = geo9.unit(name); v = len(live); G = geo9.GR
nonorth = [[s != t and G[live[s], live[t]] != 0 for t in range(v)] for s in range(v)]
rng = random.Random(seed)
def clique8():
    while True:
        order = list(range(v)); rng.shuffle(order); C = []
        for t in order:
            if all(nonorth[t][c] for c in C): C.append(t)
            if len(C) == 8: return C
C = clique8(); cpos = {t: i for i, t in enumerate(C)}
assert all(nonorth[a][b] for a in C for b in C if a != b)
top = [0]
def new():
    top[0] += 1; return top[0]
cl = []
TRUE = new(); cl.append([TRUE]); FALSE = -TRUE
def var_or_gauge(j, t):
    if t in cpos: return TRUE if cpos[t] == j else FALSE
    return new()
X = [[var_or_gauge(j, t) for t in range(v)] for j in range(K)]
Y = [[var_or_gauge(j, t) for t in range(v)] for j in range(K)]
def AND(a, b):
    if a == FALSE or b == FALSE: return FALSE
    if a == TRUE: return b
    if b == TRUE: return a
    z = new(); cl.extend([[-z, a], [-z, b], [z, -a, -b]]); return z
def XOR2(a, b):
    if a == FALSE: return b
    if b == FALSE: return a
    if a == TRUE: return -b
    if b == TRUE: return -a
    z = new(); cl.extend([[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]); return z
def XORT(lits):
    lits = list(lits)
    if not lits: return FALSE
    while len(lits) > 1:
        nxt = [XOR2(lits[i], lits[i + 1]) for i in range(0, len(lits) - 1, 2)]
        if len(lits) % 2: nxt.append(lits[-1])
        lits = nxt
    return lits[0]
def fix(lit, val):
    cl.append([lit] if val else [-lit])
nfit = 0
for s in range(v):
    for t in range(v):
        if s == t or nonorth[s][t]:
            fix(XORT(AND(Y[j][s], X[j][t]) for j in range(K)), 1 if s == t else 0); nfit += 1
z = [new() for _ in range(v)]
for j in range(K):
    fix(XORT(AND(z[t], X[j][t]) for t in range(v)), 0)
a = [XORT(z[t] for t in range(v) if (w >> t) & 1) for w in W0b]
cl.append(a)
print('unit %s K %d seed %d clique(live idx) %s -> lines %s  vars %d clauses %d fit %d solver %s' %
      (name, K, seed, C, [live[c] for c in C], top[0], len(cl), nfit, solver), flush=True)
t0 = time.time()
import threading
with Solver(name=solver, bootstrap_with=cl) as S:
    tm = threading.Timer(float(sys.argv[5]) if len(sys.argv) > 5 else 1e9, S.interrupt); tm.start()
    r = S.solve_limited(expect_interrupt=True); tm.cancel()
    print("RESULT", "UNRESOLVED (interrupted)" if r is None else "SAT (decoder with W missing W0)" if r else 'UNSAT (W contains W0 for every K-total decoder)',
          '%.1f s' % (time.time() - t0), flush=True)
    if r:
        m = set(l for l in S.get_model() if l > 0)
        val = lambda l: (l == TRUE) or (l != FALSE and l in m)
        Xs = [[int(val(X[j][t])) for t in range(v)] for j in range(K)]
        Ys = [[int(val(Y[j][t])) for t in range(v)] for j in range(K)]
        # independent check of fit
        bad = 0
        for s in range(v):
            for t in range(v):
                b = sum(Ys[j][s] * Xs[j][t] for j in range(K)) % 2
                if s == t and b != 1: bad += 1
                if s != t and nonorth[s][t] and b != 0: bad += 1
        print('fit violations', bad)
