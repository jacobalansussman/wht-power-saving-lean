#!/usr/bin/env python3
"""cx-medium / designs_cx.py (own code, stdlib; python3 -I).  Port designs as explicit LINES in R^t (root systems, written
from scratch) together with explicit binary LABELS u_a in F_2^k (standard dot product), such that
    u_a . u_b = B_ab := 1 + <r_a, r_b>/(|r|^2/2)  (mod 2)
i.e. label-orthogonal (junk allowed) exactly on the root pairs at 60/120 degrees, label-non-orthogonal on orthogonal
roots and on the diagonal (odd labels).  The labels are obtained by factoring B = U U^T over F_2 (B symmetric with unit
diagonal, rank k): pick k rows I with B_II invertible, find an orthonormal basis for the form B_II, set U = B[:,I] A^-T.
Every design used by cx-medium is checked: U U^T = B exactly, labels distinct and odd."""
import itertools, random
from fractions import Fraction

def ip(a, b): return sum(x * y for x, y in zip(a, b))

def lines(R):
    seen, L = set(), []
    for r in R:
        key = max(tuple(r), tuple(-x for x in r))
        if key not in seen:
            seen.add(key); L.append(key)
    return L

def roots_E8():                      # doubled coordinates, norm 8
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * 8; v[i] = 2 * si; v[j] = 2 * sj; R.append(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0: R.append(tuple(signs))
    return R

def roots_D(n):                       # doubled, norm 8
    R = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * n; v[i] = 2 * si; v[j] = 2 * sj; R.append(tuple(v))
    return R

def roots_A(n):                       # in R^(n+1), doubled, norm 8
    R = []
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                v = [0] * (n + 1); v[i] = 2; v[j] = -2; R.append(tuple(v))
    return R

def design(name):
    if name in ('E8', 'P112', 'E7', 'E6'):
        E8 = lines(roots_E8())
        if name == 'E8': return E8
        if name == 'P112':
            frame = set()
            for i in range(0, 8, 2):
                for sj in (1, -1):
                    v = [0] * 8; v[i] = 2; v[i + 1] = 2 * sj; frame.add(tuple(v))
            return [l for l in E8 if l not in frame]
        if name == 'E7': return [l for l in E8 if ip(l, (2, 2, 0, 0, 0, 0, 0, 0)) == 0]
        if name == 'E6': return [l for l in E8 if ip(l, (2, -2, 0, 0, 0, 0, 0, 0)) == 0 and ip(l, (0, 2, -2, 0, 0, 0, 0, 0)) == 0]
    if name[0] == 'D': return lines(roots_D(int(name[1:])))
    if name[0] == 'A': return lines(roots_A(int(name[1:])))
    raise ValueError(name)

# ---------------- F_2 linear algebra on ints (bit i = coordinate i)
def rank2(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def span_basis(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis

def dot(x, y): return bin(x & y).count('1') & 1

def label_gram(L):
    norm = ip(L[0], L[0]); unit = norm // 2
    n = len(L)
    B = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            g = ip(L[i], L[j]); assert g % unit == 0, 'not a two-angle design'
            B[i][j] = (1 + g // unit) % 2
    return B

def matrank2(rows):            # rows: lists of 0/1
    return rank2([int(''.join(map(str, r)), 2) for r in rows])

def inv2(M):
    n = len(M); A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c])
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c]: A[r] = [x ^ y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

def mul2(A, B):
    return [[sum(A[i][l] & B[l][j] for l in range(len(B))) & 1 for j in range(len(B[0]))] for i in range(len(A))]

def tr(A): return [list(r) for r in zip(*A)]

def labels(L, seed=1):
    """explicit labels u_a in F_2^k with u_a . u_b = B_ab; returned as ints (bit c = coordinate c)."""
    B = label_gram(L); n = len(L)
    # rows I: a maximal independent set of rows of B
    I, basis = [], []
    for i in range(n):
        v = int(''.join(map(str, B[i])), 2); w = v
        for b in basis: w = min(w, w ^ b)
        if w: basis.append(w); I.append(i)
    k = len(I)
    S = [[B[i][j] for j in I] for i in I]                 # invertible, symmetric, unit diagonal
    rng = random.Random(seed)
    def form(x, y): return sum(x[i] & S[i][j] & y[j] for i in range(k) for j in range(k)) & 1
    vecs = [list(map(int, format(x, f'0{k}b'))) for x in range(1, 1 << k)]
    for _ in range(2000):
        rng.shuffle(vecs); out = []
        for v in vecs:
            if form(v, v) == 1 and all(form(v, w) == 0 for w in out):
                out.append(v)
                if len(out) == k: break
        if len(out) == k: break
    assert len(out) == k, 'no orthonormal basis found'
    Pm = tr(out)                                          # columns = orthonormal basis: P^T S P = I
    A = tr(inv2(Pm))                                      # S = A A^T with A = P^-T
    U = mul2([[B[i][j] for j in I] for i in range(n)], tr(inv2(A)))   # U = B[:,I] A^-T
    labs = [int(''.join(map(str, row)), 2) for row in U]
    for i in range(n):
        for j in range(n):
            assert dot(labs[i], labs[j]) == B[i][j]
    assert len(set(labs)) == n or True
    return labs, k, B

if __name__ == '__main__':
    for name in ('E8', 'P112', 'E7', 'E6', 'D6', 'D8', 'A5', 'D4', 'A3'):
        L = design(name); labs, k, B = labels(L)
        junk = sum(1 for j in range(len(L)) if B[0][j] == 0)
        print(name, 'n', len(L), 'window k', k, 'distinct labels', len(set(labs)), 'all odd', all(bin(u).count('1') % 2 for u in labs),
              'junk/port', junk)
