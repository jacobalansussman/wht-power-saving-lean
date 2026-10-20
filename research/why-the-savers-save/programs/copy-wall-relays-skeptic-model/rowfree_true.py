"""rowfree.py (v-R1-model) -- paid variant the d-R1 pricing does not cover: RELAXED READERS.
A target S that leaves its interval (frame G not inside a_S^perp, excess exactly 2 per call: hub C.2/C.3/C.7)
can cancel any source directly, so its row of B is free.  Also allow freed COLUMNS (relayed ports T, any reader,
NO order constraint: a relaxation).  Question: rank-r F_2 decoder with |rows freed| + |cols freed| <= k ?
Gauge: k+1 pairwise disjoint r-cliques K_i; run i forbids freeing anything in K_i and sets X_Ki = I (WLOG then,
since B[K_i,K_i] = I).  By pigeonhole the union of the k+1 runs covers every choice of <= k freed rows/cols.
usage: VENV -I rowfree.py INSTANCE r k [rows|both]"""
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
import patterns
patterns.HUB_P111 = [9, 31, 37, 52, 66, 76, 78, 111, 119]   # #589 word P (e8.py order 39,52,57,71,76,78,103,109,114) mapped to hamming_e8 order
from patterns import instance
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
name, r, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
mode = sys.argv[4] if len(sys.argv) > 4 else 'rows'
NCL = int(sys.argv[5]) if len(sys.argv) > 5 else k + 1   # control runs only: fewer cliques = restriction, SAT still valid
RK = int(sys.argv[6]) if len(sys.argv) > 6 else r       # clique size (gauge X_K = [I | 0] when RK < r)
knows, G = instance(name); v = len(knows)
# disjoint cliques of pairwise non-knowing ports, greedy with restarts
import random
def clique_in(avail, size, seed):
    rnd = random.Random(seed)
    for _ in range(4000):
        cand = list(avail); rnd.shuffle(cand); cur = []
        for p in cand:
            if all(not knows[p][q] and not knows[q][p] for q in cur): cur.append(p)
            if len(cur) == size: return cur
    return None
cliques = None
for seed in range(200):
    avail = set(range(v)); cs = []
    for i in range(NCL):
        K = clique_in(avail, RK, seed * 100 + i)
        if K is None: break
        cs.append(K); avail -= set(K)
    if len(cs) == NCL: cliques = cs; break
assert cliques, 'no disjoint cliques found'
print('cliques', cliques, flush=True)
results = []
for ci, K in enumerate(cliques):
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
    fr = {S: pool.id(('fr', S)) for S in range(v)}
    fc = {T: pool.id(('fc', T)) for T in range(v)}
    Kset = set(K)
    for S in range(v):
        if S in Kset or mode == 'rows': cl.append([-fc[S]])
        if S in Kset: cl.append([-fr[S]])
    for S in range(v):
        d = Bv(S, S)
        # diagonal required unless row freed
        if d == 'F': cl.append([fr[S]])
        elif d != 'T': cl.append([fr[S], d])
        for T in range(v):
            if S == T or knows[S][T]: continue
            b = Bv(S, T)
            if b == 'F': continue
            if b == 'T': cl.append([fr[S], fc[T]])
            else: cl.append([fr[S], fc[T], -b])
    lits = [fr[S] for S in range(v)] + ([fc[T] for T in range(v)] if mode == 'both' else [])
    cl.extend(CardEnc.atmost(lits=lits, bound=k, vpool=pool, encoding=EncType.seqcounter).clauses)
    t0 = time.time()
    s = Cadical153(bootstrap_with=cl); res = s.solve()
    extra = ''
    if res:
        m = set(l for l in s.get_model() if l > 0)
        extra = 'rows freed %s cols freed %s' % ([S for S in range(v) if fr[S] in m], [T for T in range(v) if fc[T] in m])
    print('%s r=%d k=%d mode=%s clique#%d: %s  %.1fs %s' % (name, r, k, mode, ci, 'SAT' if res else 'UNSAT', time.time() - t0, extra), flush=True)
    results.append(res)
print('OVERALL', 'SAT somewhere' if any(results) else 'UNSAT in all %d runs (covers every choice)' % len(results))
