#!/usr/bin/env python3
"""v-ex2-reform: rebuild OpenAI's A0 = (J(x)J)(Gamma(x)I - I(x)Gamma) from 08-reflection.tex and recheck
ex2-reform's claim that its off-diagonal support is removal-minimal (no single pair can be dropped while some real
symmetric A with that support (diagonal free) still has kernel exactly E = ker L).  Exact arithmetic mod two primes,
own code (no numpy).  Method: basis of V_S = {A symmetric, supp in S + diag, A E = 0} once; V_{S-e} is the hyperplane
A_e = 0 inside V_S; generic rank by a random element.  Usage: python3 sparse_check.py d"""
import sys, random
d = int(sys.argv[1])
def kron(A, B):
    rb, cb = len(B), len(B[0])
    return [[A[i//rb][j//cb]*B[i%rb][j%cb] for j in range(len(A[0])*cb)] for i in range(len(A)*rb)]
def eye(n): return [[int(i==j) for j in range(n)] for i in range(n)]
def add(A, B, s=1): return [[a+s*b for a,b in zip(r1,r2)] for r1,r2 in zip(A,B)]
def mul(A, B):
    Bt = list(zip(*B)); return [[sum(a*b for a,b in zip(r,c)) for c in Bt] for r in A]
C = [[0,0],[1,0]]; Z = [[1,0],[0,-1]]; X = [[0,1],[1,0]]
D = C
for m in range(1, d): D = add(kron(C, eye(2**m)), kron(Z, D))
Gam = add(kron(X, eye(2**d)), kron(Z, D))
N = 2**(d+1); w = N*N; r = w//2
# checks on Gamma: involution, nnz
G2 = mul(Gam, Gam); assert G2 == eye(N), "Gamma not an involution"
nnzG = sum(1 for row in Gam for x in row if x)
print(f"d={d}: N={N}, w={w}; Gamma^2=I ok; nnz(Gamma)={nnzG} (paper (1+d/2)N = {(1+d/2)*N})")
J = [[1]]
for i in range(1, d+2):
    Zi = eye(2)
    for _ in range(i-1): Zi = mul(Zi, Z)
    J = kron(J, mul(Zi, X))
L = add(kron(Gam, eye(N)), kron(eye(N), Gam), -1)
S = kron(J, J)
A0 = mul(S, L)
assert all(A0[i][j] == A0[j][i] for i in range(w) for j in range(w)), "A0 not symmetric"
nnzA = sum(1 for row in A0 for x in row if x)
pairs = [(i,j) for i in range(w) for j in range(i+1,w) if A0[i][j] != 0]
diag_nz = sum(1 for i in range(w) if A0[i][i] != 0)
print(f"A0 symmetric ok; nnz(A0)={nnzA} (paper (d+2)w = {(d+2)*w}); off-diagonal pairs h={len(pairs)} (paper bound (d/2+1)w = {(d/2+1)*w}); nonzero diagonal {diag_nz}")

def nullspace_modp(rows, ncols, p):
    """rows: list of lists mod p. Returns basis of {x : rows x = 0}."""
    piv = {}  # col -> row (normalised, reduced w.r.t. earlier pivots lazily)
    order = []
    for row in rows:
        v = [x % p for x in row]
        for c in order:
            if v[c]:
                f = v[c]; pr = piv[c]
                v = [(a - f*b) % p for a,b in zip(v, pr)]
        lead = next((c for c in range(ncols) if v[c]), None)
        if lead is None: continue
        inv = pow(v[lead], p-2, p); v = [(a*inv) % p for a in v]
        # reduce existing pivots by the new one (keep RREF)
        for c in order:
            pr = piv[c]
            if pr[lead]:
                f = pr[lead]; piv[c] = [(a - f*b) % p for a,b in zip(pr, v)]
        piv[lead] = v; order.append(lead)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        x = [0]*ncols; x[fc] = 1
        for c in order: x[c] = (-piv[c][fc]) % p
        basis.append(x)
    return basis
def rank_modp(M, p):
    M = [[x % p for x in row] for row in M]; rk = 0; ncols = len(M[0])
    for c in range(ncols):
        pr = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if pr is None: continue
        M[rk], M[pr] = M[pr], M[rk]; inv = pow(M[rk][c], p-2, p)
        M[rk] = [(a*inv) % p for a in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]; M[i] = [(a - f*b) % p for a,b in zip(M[i], M[rk])]
        rk += 1
    return rk
for p in [2147483647, 1000000007]:
    rng = random.Random(12345 + p)
    # E = ker L mod p
    Kb = nullspace_modp(L, w, p)
    print(f"[p={p}] dim ker L = {len(Kb)} (need {r}); rank A0 = {rank_modp(A0, p)} (need {r}); "
          f"A0*E = 0: {all(sum(A0[i][k]*v[k] for k in range(w)) % p == 0 for v in Kb for i in range(w))}")
    cols = [(i,i) for i in range(w)] + pairs
    idx = {c:k for k,c in enumerate(cols)}
    # constraints: for each basis vector v of E and each row i: sum_j A_ij v_j = 0
    rows = []
    for v in Kb:
        for i in range(w):
            row = [0]*len(cols)
            row[idx[(i,i)]] = v[i]
            rows.append(row)
        # off-diagonal contributions
    rows2 = []
    nv = len(Kb)
    for t, v in enumerate(Kb):
        for i in range(w):
            rows[t*w + i] = rows[t*w + i]
    for k,(i,j) in enumerate(pairs):
        kk = w + k
        for t, v in enumerate(Kb):
            rows[t*w + i][kk] = (rows[t*w + i][kk] + v[j]) % p
            rows[t*w + j][kk] = (rows[t*w + j][kk] + v[i]) % p
    VS = nullspace_modp(rows, len(cols), p)
    def rand_rank(basis):
        if not basis: return 0
        coef = [rng.randrange(p) for _ in basis]
        y = [sum(c*b[k] for c,b in zip(coef, basis)) % p for k in range(len(cols))]
        A = [[0]*w for _ in range(w)]
        for k,(i,j) in enumerate(cols):
            A[i][j] = y[k]; A[j][i] = y[k]
        return rank_modp(A, p)
    full = max(rand_rank(VS) for _ in range(2))
    print(f"[p={p}] dim V_S = {len(VS)}; generic rank on full support = {full} (need {r})")
    removable = []
    for k,(i,j) in enumerate(pairs):
        kk = w + k
        # hyperplane y_kk = 0 inside span(VS)
        piv = next((b for b in VS if b[kk]), None)
        if piv is None:
            sub = VS   # coordinate already identically zero on V_S
        else:
            inv = pow(piv[kk], p-2, p)
            sub = [[(a - (b[kk]*inv) * c) % p for a,c in zip(b, piv)] for b in VS if b is not piv]
        rk = max(rand_rank(sub) for _ in range(2))
        if rk == r: removable.append((i,j))
    print(f"[p={p}] pairs whose removal keeps rank {r}: {len(removable)} of {len(pairs)} -> "
          f"{'REMOVAL-MINIMAL' if not removable else 'NOT minimal: ' + str(removable[:5])}")
