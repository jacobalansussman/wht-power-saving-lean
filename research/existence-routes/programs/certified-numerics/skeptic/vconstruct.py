# v-ex2-certnum: my own construction of a W-in-ceil(3W/2) scheme in the V(x)V stand-in (subspace-game recipe, own code),
# dumped in ex2-certnum's JSON convention (original 16-dim coordinates) for vstandin.py to verify exactly.
import sys, json
W = int(sys.argv[1]); s = -(-3 * W // 2); out = sys.argv[2]
sys.argv = sys.argv[:1]
import vstandin as V
from fractions import Fraction as Fr
Z, ONE, add, sub, mul, inv, isz, mat, mm = V.Z, V.ONE, V.add, V.sub, V.mul, V.inv, V.isz, V.mat, V.mm
idx = V.idx; Ra = V.adapted(V.R); Ka = V.adapted(V.Kr)
def rref(rows):
    M = [list(r) for r in rows]; rk = 0; piv = []
    for c in range(len(M[0])):
        p = next((i for i in range(rk, len(M)) if not isz(M[i][c])), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; iv = inv(M[rk][c]); M[rk] = [mul(iv, x) for x in M[rk]]
        for i in range(len(M)):
            if i != rk and not isz(M[i][c]):
                f = M[i][c]; M[i] = [sub(x, mul(f, y)) for x, y in zip(M[i], M[rk])]
        piv.append(c); rk += 1
    return M[:rk], piv
def blk(Ma, rl, cl): return [[Ma[r][c] for c in idx[cl]] for r in idx[rl]]
def gains(rl, cl):
    G, _ = rref(blk(Ka, rl, cl)); vs = []
    for l in range(W):
        for g in G:
            v = [Z] * (16 * W)
            for jj, c in enumerate(idx[cl]): v[16 * l + c] = g[jj]
            vs.append(v)
    return vs
GP = gains('U6', 'U4'); GQ = gains('U4', 'U6')          # P (U6-rows) needs U4-coord vectors; Q needs U6-coord vectors
RrP, pvP = rref(blk(Ra, 'U6', 'U4')); RrQ, pvQ = rref(blk(Ra, 'U4', 'U6'))
print('needs', len(GP), len(GQ), 'per-move transfer', len(RrP), len(RrQ), 'moves', s)
srcs = []
for l in range(W):
    E = mat(16, 16 * W)
    for i in range(16): E[i][16 * l + i] = ONE
    srcs.append(E)
movemix = []
for j in range(s):
    c = mat(16, 16 * W)
    for q, v in enumerate(GP[len(RrP) * j: len(RrP) * (j + 1)]): c[idx['U4'][pvP[q]]] = list(v)
    for q, v in enumerate(GQ[len(RrQ) * j: len(RrQ) * (j + 1)]): c[idx['U6'][pvQ[q]]] = list(v)
    Xs = [[[c[r][16 * l + cc] for cc in range(16)] for r in range(16)] for l in range(W)]
    movemix.append(Xs + [mat(16, 16)] * j)
    srcs.append(mm(Ra, c))
# outputs: per target t, per piece, solve row by row
def solve_rows(basis_rows, target):      # find coefficients a with sum a_i basis_rows[i] = target (exact), or None
    n = len(basis_rows); m = len(target)
    A = [[basis_rows[i][jj] for i in range(n)] + [target[jj]] for jj in range(m)]
    R_, piv = rref(A)
    if n in piv: return None
    a = [Z] * n
    for row, p in zip(R_, piv): a[p] = row[n]
    return a
outmix = []
for t in range(W):
    Xt = [mat(16, 16) for _ in srcs]
    for pc in ('U6', 'U4'):
        brs = [(i, q, srcs[i][q]) for i in range(len(srcs)) for q in idx[pc]]
        for r in idx[pc]:
            tgt = [Z] * (16 * W)
            for jj in range(16): tgt[16 * t + jj] = Ka[r][jj]
            a = solve_rows([b[2] for b in brs], tgt); assert a is not None, 'target row not reachable'
            for (i, q, _), x in zip(brs, a): Xt[i][r][q] = x
    for r in idx['L2']:
        for q in idx['L2']: Xt[t][r][q] = Ka[r][q]
    outmix.append(Xt)
# to original coordinates: M_orig = Bm M_a Bm^{-1}
Bm = [[V.BAS[c][1][r] for c in range(16)] for r in range(16)]
Bi = [[x[1] for x in V.coords([ONE if k == j else Z for k in range(16)])] for j in range(16)]   # rows: coords of e_j
Bi = [[Bi[j][r] for j in range(16)] for r in range(16)]
def orig(Ma): return mm(mm(Bm, Ma), Bi)
def js(M): return [[[str(x[0]), str(x[1])] for x in row] for row in M]
json.dump({'W': W, 's': s, 'convention': 'as ex2-certnum', 'moves': [[js(orig(M)) for M in mv] for mv in movemix],
           'outputs': [[js(orig(M)) for M in om] for om in outmix]}, open(out, 'w'))
print('wrote', out)
