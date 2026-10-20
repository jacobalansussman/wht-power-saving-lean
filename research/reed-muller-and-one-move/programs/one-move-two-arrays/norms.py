"""norms.py: total target norm sum_f W|J_f| per fibre system (pure Python, small): Singer-(2^m-1) at m = 3, 4, 5
and the lifted group at block 5.  Used only to read costs as fractions of the target."""
import sys, itertools
sys.path.insert(0, '.')
from sfield import Singer
def singer_norm(m):
    S = Singer(m, 'sd'); F = S.F; n = S.n
    seen = set(); tot = 0; nf = 0; js = []
    for Q in S.forms():
        if not any(Q): continue
        v = tuple(S.Qval(Q, x) for x in range(1, F.q))
        if v in seen: continue
        seen |= {tuple(S.Qval(Q, F.mul(F.a_exp(-e), x)) for x in range(1, F.q)) for e in range(n)}
        J = [k for k in range(n) if S.pair(Q, S.Tk(k))]
        tot += len(J); nf += 1; js.append(len(J))
    return nf, tot, sorted(set(js))
for m in (3, 4, 5):
    nf, tot, js = singer_norm(m)
    print(f"Singer m={m}: fibres {nf}, sum |J| = {tot} (per array), mean |J| {tot/nf:.2f}, |J| values {js}")

# lifted group (as in vp5.py)
S4 = Singer(4, 'sd'); F = S4.F; n = S4.n
idx5 = [(i, j) for i in range(5) for j in range(i, 5)]
def vec(x4, b5): return tuple(S4.coord[x4]) + (b5,)
def Qv(Q, v): return sum(Q[bq] & v[i] & v[j] for bq, (i, j) in enumerate(idx5)) & 1
def Tk5(k):
    T4 = S4.Tk(k); T = [[0] * 5 for _ in range(5)]
    for i in range(4):
        for j in range(4): T[i][j] = T4[i][j]
    T[4][4] = 1; return T
TK = [Tk5(k) for k in range(n)]
def pair5(Q, T): return sum(Q[bq] & T[i][j] for bq, (i, j) in enumerate(idx5)) & 1
pts = [vec(x, bb) for bb in (0, 1) for x in range(F.q) if (x, bb) != (0, 0)]
pidx = {p: i for i, p in enumerate(pts)}
perm = [pidx[vec(F.mul(2, x), bb)] for (x, bb) in [(S4.elem[p[:4]], p[4]) for p in pts]]
seen = set(); tot = 0; nf = 0; js = {}
for Q in itertools.product((0, 1), repeat=len(idx5)):
    if not any(Q): continue
    val = [Qv(Q, p) for p in pts]
    if tuple(val) in seen: continue
    cur = val
    for _ in range(n):
        seen.add(tuple(cur)); nxt = [0] * len(pts)
        for i in range(len(pts)): nxt[perm[i]] = cur[i]
        cur = nxt
    J = sum(pair5(Q, TK[k]) for k in range(n)); tot += J; nf += 1; js[J] = js.get(J, 0) + 1
print(f"lifted block 5: fibres {nf}, sum |J| = {tot} (per array), mean {tot/nf:.2f}, |J| histogram {sorted(js.items())}")
