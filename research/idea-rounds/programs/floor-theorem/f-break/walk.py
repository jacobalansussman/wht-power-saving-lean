#!/usr/bin/env python3
"""walk.py: random exact walks (mod p) testing the two claims of the transfer lemma at f = 2.
State Phi: inputs V^W -> arrays V^R, V = C^d (x) C^d.  Psi = (Phi restricted to inputs in D^W) modulo D^R.
Claim A: a cheap step (array addition with any coefficients; M (x) M on one array, M ANY monomial matrix; a diagonal
         matrix on one column; a part of an array cut by address bits scaled) never raises rank Psi.
Claim B: a block of rank k on one array raises rank Psi by at most 2^m - 2^(m-k).
Also: the plain schedule ends with rank Psi = W (2^m - 1); a shift on ONE column only can raise rank Psi (it must, to be an escape)."""
import sys, random
p = 65537; I = 256; h = pow(2, p-2, p)
def rank(rows):
    rows = [r[:] for r in rows]; rk = 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = pow(rows[rk][c], p-2, p); rows[rk] = [x*iv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]; rows[i] = [(x - f*y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def kron(A, B): return [[A[i][j]*B[k][l] % p for j in range(len(A)) for l in range(len(B))] for i in range(len(A)) for k in range(len(B))]
def mm(A, B):
    Bt = list(zip(*B)); return [[sum(a*b for a, b in zip(r, c)) % p for c in Bt] for r in A]
C = [[(1+I)*h % p, (1-I)*h % p], [(1-I)*h % p, (1+I)*h % p]]; I2 = [[1, 0], [0, 1]]
def run(m, W, R, trials, steps, seed):
    rnd = random.Random(seed); d = 2**m; n = d*d
    eye = lambda k: [[int(i == j) for j in range(k)] for i in range(k)]
    def block(k, perm):   # C on k address bits after a random relabelling of the address bits is replaced by: conjugate by a random address permutation that is LINEAR is not needed for the rank; use bits 0..k-1
        g = [[1]]
        for b in range(m): g = kron(g, C if b < k else I2)
        return g
    diag_idx = [x*d + x for x in range(d)]; off_idx = [i for i in range(n) if i not in set(diag_idx)]
    def psi_rank(Phi):  # Phi: list over arrays of n x (W*n) matrices
        rows = []
        for r in range(R):
            for i in off_idx:
                rows.append([Phi[r][i][l*n + j] for l in range(W) for j in diag_idx])
        return rank(rows)
    worstA = 0; worstB = {}; ok = True
    for t in range(trials):
        Phi = [[[int(r < W and c == r*n + i) for c in range(W*n)] for i in range(n)] for r in range(R)]
        rk = psi_rank(Phi)
        for s in range(steps):
            kind = rnd.choice(['add', 'mono', 'block', 'block', 'diag1', 'part'])
            r = rnd.randrange(R)
            if kind == 'add':
                r2 = rnd.choice([x for x in range(R) if x != r]); c = rnd.randrange(1, p); c0 = rnd.randrange(1, p)
                Phi[r] = [[(c0*a + c*b) % p for a, b in zip(ra, rb)] for ra, rb in zip(Phi[r], Phi[r2])]
                q = 0
            elif kind == 'mono':
                perm = list(range(d)); rnd.shuffle(perm); M = [[0]*d for _ in range(d)]
                for x in range(d): M[perm[x]][x] = rnd.randrange(1, p)
                Phi[r] = mm(kron(M, M), Phi[r]); q = 0
            elif kind == 'diag1':
                Dg = [[rnd.randrange(1, p) if i == j else 0 for j in range(d)] for i in range(d)]
                Phi[r] = mm(kron(Dg, eye(d)) if rnd.random() < .5 else kron(eye(d), Dg), Phi[r]); q = 0
            elif kind == 'part':     # scale the half of the array where one address bit of column 1 is 0 (and zero the rest half the time)
                bit = rnd.randrange(m); c = rnd.randrange(1, p); z = rnd.choice([0, 1])
                Phi[r] = [[(c if not ((i // d) >> bit) & 1 else z)*a % p for a in row] for i, row in enumerate(Phi[r])]; q = 0
            else:
                k = rnd.randrange(1, m+1); g = block(k, None)
                perm = list(range(d)); rnd.shuffle(perm); Pm = [[0]*d for _ in range(d)]
                for x in range(d): Pm[perm[x]][x] = 1
                Pt = [list(c) for c in zip(*Pm)]; g = mm(mm(Pm, g), Pt)      # block along an arbitrary relabelling of addresses
                Phi[r] = mm(kron(g, g), Phi[r]); q = d - d//2**k
            new = psi_rank(Phi)
            if kind in ('add', 'mono', 'diag1', 'part'):
                worstA = max(worstA, new - rk); ok &= new <= rk
            else:
                worstB[q] = max(worstB.get(q, -99), new - rk); ok &= new - rk <= q
            rk = new
    # plain schedule and a one-column shift
    Phi = [[[int(r < W and c == r*n + i) for c in range(W*n)] for i in range(n)] for r in range(R)]
    K = block(m, None)
    for r in range(W): Phi[r] = mm(kron(K, K), Phi[r])
    plain = psi_rank(Phi)
    Phi = [[[int(r < W and c == r*n + i) for c in range(W*n)] for i in range(n)] for r in range(R)]
    X = [[int(i == j ^ 1) for j in range(d)] for i in range(d)]
    Phi[0] = mm(kron(X, eye(d)), Phi[0]); sh = psi_rank(Phi)
    print(f"m={m} W={W} R={R} trials={trials} steps={steps}: all claims held={ok};  largest rise on a cheap step={worstA};"
          f"  largest rise per block by bound q: {dict(sorted(worstB.items()))};  kernel on W arrays: rank Psi={plain} (= W(2^m-1)={W*(d-1)});"
          f"  one shift of column 1 of one array: rank Psi={sh}", flush=True)
run(1, 2, 3, 300, 12, 1); run(2, 2, 3, 150, 12, 2); run(2, 1, 2, 150, 12, 3); run(3, 1, 2, 6, 8, 4)
