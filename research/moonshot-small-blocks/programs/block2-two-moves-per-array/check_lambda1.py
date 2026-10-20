"""check_lambda1.py -- v-w7-xwborder (skeptic), own code.  Lemma Lambda1 over the field F_2 (and F_3 at random), the
smallest case M = 3, W = 2, EXHAUSTIVELY over F_2:
  over T = F[eps]/(eps^2): B = B0 + eps B1 (W x M), G = G0 + eps G1 (M x W), n = n0 + eps n1 (M x M, GENERAL, not only
  strictly lower triangular), with  B G = 0  and  B n G = eps I_W.
Claim: b = rank B0, g = rank G0 satisfy W <= b + g, b <= M - W, g <= M - W (so b = g = 1 when (M, W) = (3, 2)).
Also: Proposition Lambda0's example is a solution (with triangular n); and the exact case (eps = 0 equations
B G = 0, B n G = I) has no solution at (3, 2) (Theorem 2W's linear algebra).
Method: for every (B0, G0, n0) with B0 G0 = 0 and B0 n0 G0 = 0, the remaining equations are LINEAR in (B1, G1, n1):
  B1 G0 + B0 G1 = 0,   B1 n0 G0 + B0 n1 G0 + B0 n0 G1 = I;  solvability decided by Gaussian elimination mod p."""
import itertools, sys
import numpy as np

def rank_mod(Mx, p):
    A = Mx.copy() % p; r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c] % p: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == rows: break
    return r

def solvable(Aeq, rhs, p):
    return rank_mod(Aeq, p) == rank_mod(np.c_[Aeq, rhs], p)

def linear_system(B0, G0, n0, M, W, p):
    # unknowns: B1 (W*M), G1 (M*W), n1 (M*M); equations: (N1) W*W, (T1) W*W
    nu = W * M + M * W + M * M
    rows = []; rhs = []
    def unit(k):
        v = np.zeros(nu, np.int64); v[k] = 1; return v
    iB = lambda i, j: i * M + j
    iG = lambda i, j: W * M + i * W + j
    iN = lambda i, j: W * M + M * W + i * M + j
    n0G0 = (n0 @ G0) % p; B0n0 = (B0 @ n0) % p
    for (a, c) in itertools.product(range(W), range(W)):
        # (N1)_{ac} = sum_j B1[a,j] G0[j,c] + sum_j B0[a,j] G1[j,c]
        v = np.zeros(nu, np.int64)
        for j in range(M):
            v[iB(a, j)] += G0[j, c]; v[iG(j, c)] += B0[a, j]
        rows.append(v % p); rhs.append(0)
    for (a, c) in itertools.product(range(W), range(W)):
        # (T1)_{ac} = sum_j B1[a,j] (n0 G0)[j,c] + sum_{j,k} B0[a,j] n1[j,k] G0[k,c] + sum_j (B0 n0)[a,j] G1[j,c]
        v = np.zeros(nu, np.int64)
        for j in range(M):
            v[iB(a, j)] += n0G0[j, c]; v[iG(j, c)] += B0n0[a, j]
            for k in range(M):
                v[iN(j, k)] += B0[a, j] * G0[k, c]
        rows.append(v % p); rhs.append(1 if a == c else 0)
    return np.array(rows), np.array(rhs)

def all_mats(r, c, p):
    for t in itertools.product(range(p), repeat=r * c):
        yield np.array(t, np.int64).reshape(r, c)

def exhaustive(M, W, p, tri):
    found = 0; viol = 0; seen = {}
    Bs = list(all_mats(W, M, p)); Gs = list(all_mats(M, W, p))
    if tri:
        pos = [(i, j) for i in range(M) for j in range(i)]
        ns = []
        for t in itertools.product(range(p), repeat=len(pos)):
            n = np.zeros((M, M), np.int64)
            for (i, j), v in zip(pos, t): n[i, j] = v
            ns.append(n)
    else:
        ns = list(all_mats(M, M, p))
    for B0 in Bs:
        for G0 in Gs:
            if ((B0 @ G0) % p).any(): continue
            for n0 in ns:
                if ((B0 @ n0 @ G0) % p).any(): continue
                Aeq, rhs = linear_system(B0, G0, n0, M, W, p)
                if solvable(Aeq, rhs, p):
                    found += 1
                    b, g = rank_mod(B0, p), rank_mod(G0, p)
                    seen[(b, g)] = seen.get((b, g), 0) + 1
                    if not (W <= b + g and b <= M - W and g <= M - W): viol += 1
    return found, viol, seen

if __name__ == "__main__":
    # (1) exact case control at (3,2) over F_2: B G = 0, B n G = I impossible (general n)
    p = 2; M, W = 3, 2; exact = 0
    for B in all_mats(W, M, p):
        for G in all_mats(M, W, p):
            if ((B @ G) % p).any(): continue
            for n in all_mats(M, M, p):
                if (((B @ n @ G) % p) == np.eye(W, dtype=np.int64)).all(): exact += 1
    print("exact case (3,2) over F_2, general n: solutions =", exact, "(Theorem 2W linear algebra: expect 0)", flush=True)
    for tri in (True, False):
        found, viol, seen = exhaustive(3, 2, 2, tri)
        print("order-1 (3,2) over F_2, n %s: (B0,G0,n0) with a solution: %d; Lambda1 violations: %d; (b,g) counts %s"
              % ("strictly lower triangular" if tri else "general", found, viol, seen), flush=True)
    # (2) (M, W) = (4, 3) over F_2 at order 1 must be impossible (M >= 3W/2 = 4.5): exhaustive over B0, G0 in
    #     normal form G0 = [I_g 0; 0 0] (left GL_M and right GL_W act on the equations), every n0 (strictly lower
    #     triangular after the change of basis is NOT preserved, so n general: 2^16) -- too big; do g-normal form
    #     with B0 general (2^12) and n0 general (2^16) only when B0 G0 = 0 and the pair passes the rank pre-test.
    M, W, p = 4, 3, 2
    cnt = 0; solv = 0
    pairs = []
    for g in range(0, W + 1):
        G0 = np.zeros((M, W), np.int64)
        for i in range(g): G0[i, i] = 1
        for B0 in all_mats(W, M, p):
            if ((B0 @ G0) % p).any(): continue
            b = rank_mod(B0, p)
            pairs.append((B0, G0, b, g))
    print("(4,3): %d (B0, G0) pairs with G0 in normal form and B0 G0 = 0" % len(pairs), flush=True)
    # necessary condition from the proof (b + g >= W and b, g <= 1): no pair can pass -> report how many pass
    ok = [(b, g) for (_, _, b, g) in pairs if W <= b + g and b <= M - W and g <= M - W]
    print("(4,3): pairs meeting Lambda1's necessary rank conditions:", len(ok), "(expect 0)", flush=True)
    # brute force on a random subset of n0 for the pairs with the largest b + g, to confirm unsolvability directly
    rng = np.random.default_rng(3); tested = 0; hits = 0
    big = sorted(pairs, key=lambda t: -(t[2] + t[3]))[:60]
    for (B0, G0, b, g) in big:
        for t in range(300):
            n0 = rng.integers(0, 2, (M, M))
            if ((B0 @ n0 @ G0) % p).any(): continue
            Aeq, rhs = linear_system(B0, G0, n0, M, W, p); tested += 1
            if solvable(Aeq, rhs, p): hits += 1
    print("(4,3): direct solvability tests on %d (B0,G0,n0) triples with the largest b+g: %d solvable (expect 0)"
          % (tested, hits), flush=True)
