# construct.py W s -- EXACT construction (sympy, Gaussian rationals) of a W-in-s scheme in the V(x)V stand-in, s = ceil(3W/2),
# from the structure found in decomp.py / sym2info.py, and its exact verification in the ORIGINAL 16-dim coordinates.
#  Lemma 1 (checked exactly here): A = span rho(B) = End(L2) (+) End(U6) (+) End(U4), where U6 = off-diagonal symmetric
#    tensors, U4 = diagonal tensors e_a e_a, L2 = antisymmetric tensors.  (Every rho(b) is block diagonal in this basis,
#    and dim A = 88 = 36 + 36 + 16, so every block-diagonal matrix lies in A.)
#  Lemma 2 (the subspace game): with P = sum of row spaces of the U6-rows of the sources and Q the same for the U4-rows, a
#    move adds to P any 2 vectors of the old Q (rank R_{U6,U4} = 2) and to Q any 2 vectors of the old P (rank R_{U4,U6} = 2);
#    the outputs need P >= (U6 + U6 K) (x) C^W and Q >= (U4 + U4 K) (x) C^W.
#  Construction: the needed gains (U6+U6K) cap U4 (3-dim per array) and (U4+U4K) cap U6 (3-dim per array) are already in
#    the initial Q and P, so s = ceil(3W/2) moves on the INPUTS ALONE transfer them two at a time; the Lambda^2 part of
#    every target already lies in A (output t = K on input t there).  Outputs: exact linear solves.
import sys, itertools, time
import sympy as sp
W, s = int(sys.argv[1]), int(sys.argv[2]); t0 = time.time()
I = sp.I; n = 4
def perm(f):
    M = sp.zeros(n, n)
    for x in range(n): M[f(x), x] = 1
    return M
X0 = perm(lambda x: x ^ 1); X1 = perm(lambda x: x ^ 2)
CN01 = perm(lambda x: x ^ ((x & 1) << 1)); CN10 = perm(lambda x: x ^ ((x >> 1) & 1))
S0 = sp.diag(*[I if (x & 1) else 1 for x in range(n)]); S1 = sp.diag(*[I if (x & 2) else 1 for x in range(n)])
gens = [X0, X1, CN01, CN10, S0, S1, I * sp.eye(n)]
ca, cb = (1 + I) / 2, (1 - I) / 2
D0 = ca * sp.eye(n) + cb * X0; D1 = ca * sp.eye(n) + cb * X1; K = sp.expand(D1 * D0)
def kron(A, B): return sp.Matrix(sp.BlockMatrix([[A[i, j] * B for j in range(A.cols)] for i in range(A.rows)]).as_explicit())
rho = lambda g: sp.expand(kron(g, g))
R = rho(D0); Kr = rho(K); d = 16
# adapted basis of V(x)V: U6 (a<b symmetric), U4 (diagonal), L2 (a<b antisymmetric) -- columns of Bm
cols = []
for a in range(4):
    for b in range(a + 1, 4):
        v = sp.zeros(16, 1); v[a * 4 + b] = 1; v[b * 4 + a] = 1; cols.append(v)
for a in range(4):
    v = sp.zeros(16, 1); v[a * 4 + a] = 1; cols.append(v)
for a in range(4):
    for b in range(a + 1, 4):
        v = sp.zeros(16, 1); v[a * 4 + b] = 1; v[b * 4 + a] = -1; cols.append(v)
Bm = sp.Matrix.hstack(*cols); Bi = Bm.inv()
blocks = [range(0, 6), range(6, 10), range(10, 16)]
def blockdiag_ok(M):
    Mp = Bi * M * Bm
    return all(Mp[i, j] == 0 for i in range(16) for j in range(16) if not any(i in b and j in b for b in blocks))
# Lemma 1, exact: the generators are block diagonal (then so is the whole group B and A), and dim A = 88
assert all(blockdiag_ok(rho(g)) for g in gens), 'a generator is not block diagonal'
# dim A: BFS the group (as exact matrices with entries in Z[i]) and take the rank of the span of rho(B) exactly
seen = {tuple(sp.eye(n)): sp.eye(n)}; frontier = [sp.eye(n)]
while frontier:
    nf = []
    for g in frontier:
        for h in gens:
            M = h * g; k_ = tuple(M)
            if k_ not in seen: seen[k_] = M; nf.append(M)
    frontier = nf
print('|B| =', len(seen), '(%.0fs)' % (time.time() - t0), flush=True)
# rank of span{rho(b)} over Q(i): realify (re, im parts) and use exact integer rank via flint if present, else sympy
vecs = {}
for g in seen.values():
    v = tuple(sp.expand(x) for x in rho(g)); vecs[v] = 1
rows = [[int(sp.re(x)) for x in v] + [int(sp.im(x)) for x in v] for v in vecs] + \
       [[-int(sp.im(x)) for x in v] + [int(sp.re(x)) for x in v] for v in vecs]
try:
    from flint import fmpz_mat
    rA = fmpz_mat(rows).rank() // 2
except ImportError:
    rA = sp.Matrix(rows).rank() // 2
print('distinct rho(b): %d; dim_C A = %d (exact rank; 36 + 36 + 16 = 88 needed for Lemma 1)' % (len(vecs), rA), flush=True)
assert rA == 88
Rp = Bi * R * Bm; Kp = Bi * Kr * Bm
def symalt_ok(M):
    Mp = Bi * M * Bm
    return all(Mp[i, j] == 0 for i in range(16) for j in range(16) if (i < 10) != (j < 10))
assert symalt_ok(R) and symalt_ok(Kr)             # R and K preserve Sym^2 and Lambda^2 (they mix U6 and U4)
U6, U4, L2 = list(range(0, 6)), list(range(6, 10)), list(range(10, 16))
# Lemma 2 data
def sub(M, r, c): return M.extract(r, c)
print('ranks: R[U6,U4] =', sub(Rp, U6, U4).rank(), ' R[U4,U6] =', sub(Rp, U4, U6).rank(), flush=True)
def rowspace(M): return [M.row(i) for i in range(M.rows)]
def basis_of(vs):
    M = sp.Matrix.vstack(*vs); return [M.row(i) for i in range(M.rows)] if M.rank() == M.rows else [r for r in M.rref()[0].tolist() if any(x != 0 for x in r)]
Ksym = sub(Kp, U6 + U4, U6 + U4); Rsym = sub(Rp, U6 + U4, U6 + U4)
e = sp.eye(10)
def cap_with_coords(span_rows, coords):
    # basis of {x in span(span_rows) : x supported on coords}: solve for combinations vanishing off coords
    M = sp.Matrix.vstack(*span_rows); off = [j for j in range(10) if j not in coords]
    N = sp.Matrix(M.extract(list(range(M.rows)), off)).T.nullspace()
    return [sp.expand((v.T * M)) for v in N]
HP = [e.row(i) for i in range(6)] + [Ksym.row(i) for i in range(6)]       # U6 + U6 K  (row vectors)
LQ = [e.row(i) for i in range(6, 10)] + [Ksym.row(i) for i in range(6, 10)]  # U4 + U4 K
gP = cap_with_coords(HP, list(range(6, 10))); gQ = cap_with_coords(LQ, list(range(0, 6)))
gP = [g for g in sp.Matrix.vstack(*gP).rref()[0].tolist() if any(x != 0 for x in g)]
gQ = [g for g in sp.Matrix.vstack(*gQ).rref()[0].tolist() if any(x != 0 for x in g)]
print('dim (U6+U6K) =', sp.Matrix.vstack(*HP).rank(), ' needed P-gain per array dim', len(gP),
      ';  dim (U4+U4K) =', sp.Matrix.vstack(*LQ).rank(), ' needed Q-gain per array dim', len(gQ), flush=True)
need = -(-len(gP) * W // 2), -(-len(gQ) * W // 2)
print('lower bound from Lemma 2: s >= %d (P) and %d (Q); constructing with s = %d' % (need[0], need[1], s), flush=True)
# needed gain vectors in C^{10W}: g (x) e_l
def place(g, l):
    v = [0] * (10 * W)
    for j in range(10): v[l * 10 + j] = g[j]
    return v
GP = [place(g, l) for l in range(W) for g in gP]; GQ = [place(g, l) for l in range(W) for g in gQ]
# pick independent columns of R[U6,U4] and R[U4,U6]
RUU4 = sub(Rsym, list(range(6)), list(range(6, 10))); RU4U = sub(Rsym, list(range(6, 10)), list(range(6)))
pc = RUU4.rref()[1][:2]; qc = RU4U.rref()[1][:2]
# sources in Sym^2 coordinates (10 x 10W); inputs E_l
E = []
for l in range(W):
    M = sp.zeros(10, 10 * W); M[:, l * 10:(l + 1) * 10] = sp.eye(10); E.append(M)
srcs = list(E); mix = []      # mix[j] = list over sources of (X6, Y4)
for j in range(s):
    c = sp.zeros(10, 10 * W)
    ys = GP[2 * j:2 * j + 2]; zs = GQ[2 * j:2 * j + 2]
    for k_, y in enumerate(ys):           # c_U4 row (6 + qcol) = y   (y lies in the initial Q = U4 (x) C^W)
        for jj in range(10 * W): c[6 + pc[k_], jj] = y[jj]
    for k_, z in enumerate(zs):           # c_U6 row qcol = z         (z lies in the initial P = U6 (x) C^W)
        for jj in range(10 * W): c[qc[k_], jj] = z[jj]
    # write c = sum_l diag(X_l, Y_l) E_l: X_l = c[U6, block l][:, U6 coords], etc. (c is supported correctly)
    Xs = []
    for l in range(W):
        blk = c[:, l * 10:(l + 1) * 10]
        assert all(blk[i, jj] == 0 for i in range(6) for jj in range(6, 10)) and all(blk[i, jj] == 0 for i in range(6, 10) for jj in range(6))
        Xs.append((blk.extract(list(range(6)), list(range(6))), blk.extract(list(range(6, 10)), list(range(6, 10)))))
    Xs += [(sp.zeros(6, 6), sp.zeros(4, 4))] * (len(srcs) - W)
    mix.append(Xs); srcs.append(sp.expand(Rsym * c))
# outputs: solve exactly for (X, Y) on all sources with sum_i diag(X_i, Y_i) s_i = K in block t
outmix = []
for t in range(W):
    Tt = sp.zeros(10, 10 * W); Tt[:, t * 10:(t + 1) * 10] = Ksym
    Xt = [sp.zeros(6, 6) for _ in srcs]; Yt = [sp.zeros(4, 4) for _ in srcs]
    for part, rows_, size in [('U6', list(range(6)), 6), ('U4', list(range(6, 10)), 4)]:
        Ms = sp.Matrix.vstack(*[sp.Matrix(si.extract(rows_, list(range(10 * W)))) for si in srcs])   # (size*#src) x 10W
        for r_ in range(size):
            target = Tt.row(rows_[r_])
            sol, params = sp.Matrix(Ms.T).gauss_jordan_solve(target.T)
            sol = sol.subs({p_: 0 for p_ in params})
            for i in range(len(srcs)):
                for q in range(size):
                    if part == 'U6': Xt[i][r_, q] = sol[i * size + q]
                    else: Yt[i][r_, q] = sol[i * size + q]
    outmix.append(list(zip(Xt, Yt)))
print('Sym^2 scheme built (%.0fs)' % (time.time() - t0), flush=True)
# ---- exact verification in the ORIGINAL 16-dim coordinates: lift every mixing matrix to V(x)V (block diagonal in the
# adapted basis, with Lambda^2 part 0 for moves and K|L2 on input t for output t), then run the program exactly.
KL2 = sub(Kp, L2, L2)
def lift(X6, Y4, L):
    Mp = sp.zeros(16, 16)
    Mp[0:6, 0:6] = X6; Mp[6:10, 6:10] = Y4; Mp[10:16, 10:16] = L
    return sp.expand(Bm * Mp * Bi)
E16 = []
for l in range(W):
    M = sp.zeros(16, 16 * W); M[:, l * 16:(l + 1) * 16] = sp.eye(16); E16.append(M)
src16 = list(E16)
for j in range(s):
    c = sp.zeros(16, 16 * W)
    for i, (X6, Y4) in enumerate(mix[j]): c += lift(X6, Y4, sp.zeros(6, 6)) * src16[i]
    src16.append(sp.expand(R * c))
ok = True; nnz = 0
for t in range(W):
    o = sp.zeros(16, 16 * W)
    for i, (X6, Y4) in enumerate(outmix[t]):
        L = KL2 if i == t else sp.zeros(6, 6)
        Mi = lift(X6, Y4, L); nnz += sum(1 for x in Mi if x != 0)
        assert blockdiag_ok(Mi)
        o += Mi * src16[i]
    target = sp.zeros(16, 16 * W); target[:, t * 16:(t + 1) * 16] = Kr
    diff = sp.expand(o - target); bad = sum(1 for x in diff if x != 0)
    print('output %d: %d of %d entries differ from rho(D1 D0) in block %d (exact)' % (t, bad, len(diff), t), flush=True)
    ok &= bad == 0
print('EXACT VERDICT (%.0fs): %s' % (time.time() - t0, ('a %d-in-%d scheme EXISTS in the V(x)V stand-in: every mixing matrix is block '
      'diagonal in the adapted basis, hence in A by Lemma 1, and the program output equals the target exactly' % (W, s)) if ok else 'FAILED'), flush=True)
# dump the exact scheme in the ORIGINAL 16-dim coordinates for an independent check: every mixing matrix as a 16 x 16
# matrix of Gaussian rationals [["re", "im"], ...], moves[j][i] = coefficient of source i in move j (sources: inputs
# 0..W-1, then move outputs W..W+s-1), outputs[t][i] likewise; R = rho(D0), target rho(D1 D0) in block t.
import json
def mjson(M): return [[[str(sp.re(x)), str(sp.im(x))] for x in M.row(i)] for i in range(M.rows)]
moves_j = [[mjson(lift(X6, Y4, sp.zeros(6, 6))) for (X6, Y4) in mix[j]] for j in range(s)]
outs_j = [[mjson(lift(X6, Y4, KL2 if i == t else sp.zeros(6, 6))) for i, (X6, Y4) in enumerate(outmix[t])] for t in range(W)]
json.dump({'W': W, 's': s, 'convention': 'sources are 16 x 16W; move j: src_{W+j} = R @ sum_i moves[j][i] @ src_i; '
           'output t: sum_i outputs[t][i] @ src_i must equal rho(D1 D0) in block t; rho(g) = kron(g, g); '
           'D0 = (1+i)/2 I + (1-i)/2 X0 with X0 = x -> x^1 on C^4 (basis index x = 2*bit1 + bit0)',
           'moves': moves_j, 'outputs': outs_j}, open('scheme_%din%d.json' % (W, s), 'w'))
print('wrote scheme_%din%d.json' % (W, s), flush=True)
