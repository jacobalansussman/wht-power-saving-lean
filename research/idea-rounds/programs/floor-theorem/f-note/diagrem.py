# First family of the proposed experiment (f-note.md section 4), length N = 16:
# remainders D = 4 * (diagonal of signs), all 2^16 of them.
# (a) for which sign patterns does H - D have rank <= N/2 = 8 ?  (rank modulo a prime, then exact for the hits)
# (b) for those: is the column space of H - D spanned by vectors with at most 3 non-zero entries ?
#     A vector supported on S lies in the column space iff L[:, S] v = 0, L = a basis of the left kernel.
import sys
from fractions import Fraction
from itertools import combinations
P = 998244353
N = 16
H = [[(-1) ** bin(x & y).count("1") for y in range(N)] for x in range(N)]
def rank_mod(M):
    M = [r[:] for r in M]; r = 0
    for c in range(len(M[0])):
        p = next((i for i in range(r, len(M)) if M[i][c] % P), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = pow(M[r][c], P - 2, P)
        for i in range(r + 1, len(M)):
            if M[i][c] % P:
                f = M[i][c] * inv % P
                M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r])]
        r += 1
    return r
def rref_q(M):
    M = [[Fraction(v) for v in r] for r in M]; r = 0; piv = []
    for c in range(len(M[0])):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        d = M[r][c]; M[r] = [a / d for a in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    return M[:r], piv
def left_kernel(R):   # rows w with w R = 0, exact
    T = [[R[x][y] for x in range(N)] for y in range(N)]      # transpose: T w = 0
    rows, piv = rref_q(T)
    free = [c for c in range(N) if c not in piv]
    out = []
    for fcol in free:
        w = [Fraction(0)] * N; w[fcol] = Fraction(1)
        for r, pc in zip(rows, piv): w[pc] = -r[fcol]
        out.append(w)
    return out
hist = {}; hits = []
for bits in range(1 << N):
    R = [[H[x][y] - (4 * (1 if (bits >> x) & 1 else -1) if x == y else 0) for y in range(N)] for x in range(N)]
    rk = rank_mod(R)
    hist[rk] = hist.get(rk, 0) + 1
    if rk <= 8: hits.append(bits)
print("rank of H - 4*signs modulo", P, ": histogram", sorted(hist.items()))
print("sign patterns with rank <= 8:", len(hits))
done = 0; good = 0
for bits in hits[:400]:
    R = [[H[x][y] - (4 * (1 if (bits >> x) & 1 else -1) if x == y else 0) for y in range(N)] for x in range(N)]
    rows, piv = rref_q(R); rk = len(piv)
    L = left_kernel(R)
    sparse = []
    for k in (1, 2, 3):
        for S in combinations(range(N), k):
            cols = [[w[y] for y in S] for w in L]
            rr, pv = rref_q(cols)
            if len(pv) < k: sparse.append(S)
    # dimension of the span of the sparse vectors: build one kernel vector per dependent support
    vecs = []
    for S in sparse:
        cols = [[w[y] for y in S] for w in L]
        rr, pv = rref_q(cols)
        free = [c for c in range(len(S)) if c not in pv]
        v = [Fraction(0)] * N
        fc = free[0]; v[S[fc]] = Fraction(1)
        for r, pc in zip(rr, pv): v[S[pc]] = -r[fc]
        vecs.append(v)
    dim = len(rref_q(vecs)[1]) if vecs else 0
    done += 1
    if dim >= rk and rk <= 8: good += 1
    if done <= 6 or (dim >= rk and rk <= 8):
        print(" signs %04x: exact rank %d, supports of size <= 3 carrying a column-space vector: %d, their span: %d" % (bits, rk, len(sparse), dim))
print("patterns examined exactly:", done, "; with a basis of vectors of at most 3 non-zeros:", good)
