"""g2 (g-codes): search for a SYMMETRIC rank-8 decoder B = X X^T over F_p (p odd) on the bit E8 pattern:
x_S in F_p^8 supported on Z(S) (gauge K = I), x_S . x_T = 0 for every non-orthogonal pair, x_S . x_S != 0.
Backtracking with forward checking (domains = projective vectors on Z(S)).  usage: g2_rank8_fp.py p [nodes]"""
import sys, itertools, time
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
from copyfloor import e8_lines
p = int(sys.argv[1]); LIM = int(sys.argv[2]) if len(sys.argv) > 2 else 2_000_000
L, G = e8_lines(); n = 120
nonorth = [[S != T and G[S][T] != 0 for T in range(n)] for S in range(n)]
K = (0, 8, 9, 16, 17, 24, 25, 32)
idx = {k: i for i, k in enumerate(K)}
def proj(Zi):
    out = []
    for t in itertools.product(range(p), repeat=len(Zi)):
        if any(t) and t[next(i for i in range(len(t)) if t[i])] == 1:
            v = [0] * 8
            for c, i in zip(t, Zi): v[i] = c
            if sum(a * a for a in v) % p: out.append(tuple(v))
    return out
dom = {}
for S in range(n):
    if S in K:
        v = [0] * 8; v[idx[S]] = 1; dom[S] = [tuple(v)]
    else:
        dom[S] = proj([idx[k] for k in K if G[S][k] == 0])
order = sorted(range(n), key=lambda S: (len(dom[S]), S))
assign = {}; nodes = [0]; t0 = time.time()
def dot(a, b): return sum(x * y for x, y in zip(a, b)) % p
def rec(i, doms):
    nodes[0] += 1
    if nodes[0] > LIM: raise TimeoutError
    if i == n: return True
    # pick unassigned with smallest domain
    S = min((T for T in range(n) if T not in assign), key=lambda T: len(doms[T]))
    for x in doms[S]:
        nd = {}
        ok = True
        for T in range(n):
            if T in assign or T == S: continue
            if nonorth[S][T]:
                d = [y for y in doms[T] if dot(x, y) == 0]
                if not d: ok = False; break
                nd[T] = d
            else:
                nd[T] = doms[T]
        if not ok: continue
        assign[S] = x
        if rec(i + 1, nd): return True
        del assign[S]
    return False
try:
    r = rec(0, dom)
    print('p=%d symmetric rank-8 on bit E8 pattern: %s  nodes %d  %.1fs' % (p, 'FOUND' if r else 'UNSAT (exhaustive, gauge K fixed)', nodes[0], time.time() - t0))
    if r:
        B = [[dot(assign[S], assign[T]) for T in range(n)] for S in range(n)]
        bad = [(S, T) for S in range(n) for T in range(n) if (S == T and B[S][S] == 0) or (nonorth[S][T] and B[S][T])]
        print('check bad entries:', len(bad))
except TimeoutError:
    print('p=%d: node limit %d reached, undecided (%.1fs)' % (p, LIM, time.time() - t0))
