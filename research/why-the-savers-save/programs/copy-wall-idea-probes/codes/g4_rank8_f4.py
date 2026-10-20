"""g4 (g-codes, F_4 = vector-linear F_2 with t = 2 restricted to F_4-linear): ASYMMETRIC rank-8 decoder B = X Y^T over F_p on the bit E8 pattern (general gauge: B[K,K] = I
after row scaling, X_K = Y_K = I).  x_S, y_S in F_p^8 supported on Z(S); x_S.y_T = 0 = x_T.y_S for non-orthogonal
S != T; x_S.y_S != 0.  Backtracking with forward checking; exhaustive => no rank-8 decoder over F_p at all
(every rank-8 fitting B has an 8-clique block = diagonal, and the gauge is without loss of generality).
usage: g3_rank8_fp_asym.py p [node limit] [clique index]"""
import sys, itertools, time
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
from copyfloor import e8_lines
p = 4; LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
ADD = [[a ^ b for b in range(4)] for a in range(4)]
MUL = [[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]
L, G = e8_lines(); n = 120
nonorth = [[S != T and G[S][T] != 0 for T in range(n)] for S in range(n)]
K = (0, 8, 9, 16, 17, 24, 25, 32)
for a, b in itertools.combinations(K, 2): assert nonorth[a][b]
idx = {k: i for i, k in enumerate(K)}
def proj(Zi):
    out = []
    for t in itertools.product(range(p), repeat=len(Zi)):
        if any(t) and t[next(i for i in range(len(t)) if t[i])] == 1:
            v = [0] * 8
            for c, i in zip(t, Zi): v[i] = c
            out.append(tuple(v))
    return out
def dot(a, b):
    s = 0
    for x, y in zip(a, b): s = ADD[s][MUL[x][y]]
    return s
DX, DY = {}, {}
for S in range(n):
    if S in K:
        v = [0] * 8; v[idx[S]] = 1; DX[S] = [tuple(v)]; DY[S] = [tuple(v)]
    else:
        DX[S] = proj([idx[k] for k in K if G[S][k] == 0]); DY[S] = list(DX[S])
assign = {}; nodes = [0]; t0 = time.time()
def viable(dx, dy):
    for x in dx:
        for y in dy:
            if dot(x, y): return True
    return False
def rec(DX, DY):
    nodes[0] += 1
    if nodes[0] > LIM: raise TimeoutError
    free = [T for T in range(n) if T not in assign]
    if not free: return True
    S = min(free, key=lambda T: len(DX[T]) * len(DY[T]))
    for x in DX[S]:
        for y in DY[S]:
            if not dot(x, y): continue
            nX, nY = dict(DX), dict(DY); ok = True
            for T in free:
                if T == S or not nonorth[S][T]: continue
                dx = [u for u in DX[T] if dot(u, y) == 0]
                dy = [w for w in DY[T] if dot(x, w) == 0]
                if not dx or not dy or not viable(dx, dy): ok = False; break
                nX[T], nY[T] = dx, dy
            if not ok: continue
            assign[S] = (x, y)
            if rec(nX, nY): return True
            del assign[S]
    return False
try:
    r = rec(DX, DY)
    print('p=%d asymmetric rank-8 on bit E8 pattern: %s  nodes %d  %.1fs' % (p, 'FOUND' if r else 'UNSAT (exhaustive)', nodes[0], time.time() - t0))
    if r:
        bad = 0
        for S in range(n):
            for T in range(n):
                v = dot(assign[S][0], assign[T][1])
                if (S == T and v == 0) or (nonorth[S][T] and v): bad += 1
        print('verified, bad entries:', bad)
        import json
        json.dump({str(S): assign[S] for S in assign}, open('out/g3_rank8_p%d.json' % p, 'w'))
except TimeoutError:
    print('p=%d: node limit %d reached, undecided (%.1fs)' % (p, LIM, time.time() - t0))
