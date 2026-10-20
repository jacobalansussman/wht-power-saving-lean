"""g5 (g-codes): rank-8 decoder over F_p on the bit P111 pattern (E8 minus the 9 punctured lines of copyfloor's
calibration).  Same gauge and search as g3.  usage: g5_p111.py p [node limit]"""
import sys, itertools, time
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
from copyfloor import e8_lines
p = int(sys.argv[1]); LIM = int(sys.argv[2]) if len(sys.argv) > 2 else 300000
L, G0 = e8_lines()
P = [8, 9, 10, 12, 16, 17, 24, 32, 64]
live = [i for i in range(120) if i not in P]; n = len(live)
G = [[G0[a][b] for b in live] for a in live]
nonorth = [[S != T and G[S][T] != 0 for T in range(n)] for S in range(n)]
adj = [sum(1 << T for T in range(n) if nonorth[S][T]) for S in range(n)]
def clique8():
    def rec(cur, cand):
        if len(cur) == 8: return cur
        c = cand
        while c:
            low = c & -c; i = low.bit_length() - 1; c ^= low
            r = rec(cur + [i], cand & adj[i] & ~((low << 1) - 1))
            if r: return r
        return None
    return rec([], (1 << n) - 1)
K = tuple(clique8()); print('live lines', n, 'clique', K)
idx = {k: i for i, k in enumerate(K)}
def proj(Zi):
    out = []
    for t in itertools.product(range(p), repeat=len(Zi)):
        if any(t) and t[next(i for i in range(len(t)) if t[i])] == 1:
            v = [0] * 8
            for c, i in zip(t, Zi): v[i] = c
            out.append(tuple(v))
    return out
def dot(a, b): return sum(x * y for x, y in zip(a, b)) % p
DX, DY = {}, {}
for S in range(n):
    if S in K:
        v = [0] * 8; v[idx[S]] = 1; DX[S] = [tuple(v)]; DY[S] = [tuple(v)]
    else:
        DX[S] = proj([idx[k] for k in K if G[S][k] == 0]); DY[S] = list(DX[S])
assign = {}; nodes = [0]; t0 = time.time()
def viable(dx, dy):
    return any(dot(x, y) for x in dx for y in dy)
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
    print('p=%d rank-8 on bit P111 pattern: %s  nodes %d  %.1fs' % (p, 'FOUND' if r else 'UNSAT (exhaustive)', nodes[0], time.time() - t0))
except TimeoutError:
    print('p=%d: node limit reached, undecided (%.1fs)' % (p, time.time() - t0))
