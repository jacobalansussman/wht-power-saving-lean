# Exact check of the note writer's added statement (f-note.md section 3):
# every non-zero vector v with H v = -sqrt(N) v has at least sqrt(N) non-zero entries,
# and rank(H - sqrt(N) I) = N/2.  H = plain +-1 Walsh-Hadamard matrix, N = 4, 16.
# A vector supported on S lies in that eigenspace iff the columns S of (H + sqrt(N) I) are dependent.
from fractions import Fraction
from itertools import combinations
def had(N):
    return [[(-1) ** bin(x & y).count("1") for y in range(N)] for x in range(N)]
def rank(rows):
    M = [[Fraction(v) for v in r] for r in rows]
    r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r
for N, s in ((4, 2), (16, 4)):
    H = had(N)
    minus = [[H[x][y] - (s if x == y else 0) for y in range(N)] for x in range(N)]
    plus = [[H[x][y] + (s if x == y else 0) for y in range(N)] for x in range(N)]
    print("N =", N, " rank(H - sqrt(N) I) =", rank(minus), " rank(H + sqrt(N) I) =", rank(plus))
    for k in range(1, s + 1):
        dep = 0
        tot = 0
        for S in combinations(range(N), k):
            tot += 1
            cols = [[plus[x][y] for y in S] for x in range(N)]
            if rank(cols) < k:
                dep += 1
        print("  supports of size", k, ":", tot, "sets,", dep, "carry a (-sqrt N)-eigenvector")
# control: the same test must find eigenvectors at size sqrt(N) (it does if the last line shows dep > 0)
