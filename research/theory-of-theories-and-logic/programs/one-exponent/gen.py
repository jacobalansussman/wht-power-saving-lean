# b56-methods gen.py: the dual polar graph of Sp(2m,2) in the self-dual-code model, and a starting lift.
# Model: n = 2m+2 points, point INF = n-1.  Symplectic label space = even subsets of n points modulo the all-ones
# word, form |A cap B| mod 2 (nondegenerate, dim 2m).  Lagrangians <-> self-dual binary codes C of length n
# (dim m+1, contain all-ones).  Punctured model: L = {c in C : INF not in c}, an m-dim space of even subsets of
# the first n-1 = 2m+1 points (this is the w-bounds "seven-point model" at m = 3).
# Lift: m integer rows of length 2m+1 whose reductions mod 2 form a basis of L.  Pair frames (codes spanned by
# n/2 disjoint pairs) get rows e_a - e_b for their pairs avoiding INF: delta = distance exactly (perfect).
# Other frames start from their 0/1 echelon basis.
# Usage: python3 -I gen.py m out.txt   -> writes: nframes m d deg ; per frame: fixed flag, m*d ints, deg neighbours
import sys, itertools
m = int(sys.argv[1]); n = 2 * m + 2; INF = n - 1; ONE = (1 << n) - 1; d = n - 1
def rref(vecs):
    basis = []  # list of (pivot, vec)
    for v in vecs:
        for p, b in basis:
            if v >> p & 1: v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [(q, (b ^ v) if (b >> p & 1) else b) for q, b in basis]
            basis.append((p, v))
    basis.sort(reverse=True)
    return tuple(b for p, b in basis)
def span(basis):
    out = [0]
    for b in basis: out += [x ^ b for x in out]
    return out
def dot(a, b): return bin(a & b).count('1') & 1
def perp(basis):
    # all vectors of F2^n orthogonal to basis: solve via brute force on a complement (n <= 16 fine)
    B = rref(basis); piv = [b.bit_length() - 1 for b in B]
    free = [i for i in range(n) if i not in piv]
    sol = []
    for f in free:
        v = 1 << f
        # choose pivot bits so that dot with each basis vector is 0; B is in RREF with distinct pivots
        for b, p in zip(B, piv):
            if dot(v, b): v ^= 1 << p
        sol.append(v)
    # verify
    for v in sol:
        for b in B: assert dot(v, b) == 0
    return rref(sol)
start = rref([ONE] + [(1 << (2 * i)) | (1 << (2 * i + 1)) for i in range(m + 1)])
assert len(start) == m + 1
idx = {start: 0}; codes = [start]; nbrs = []
q = 0
while q < len(codes):
    C = codes[q]; q += 1
    words = span(C)
    # hyperplanes containing ONE: kernels of nonzero functionals on C/<ONE>; get them as {c : c.h = 0} for h in a
    # complement... simpler: enumerate (m)-dim subspaces C0 with ONE in C0 via functionals on a basis.
    Cb = list(C)
    # basis with ONE first
    B = [ONE];
    for b in Cb:
        if len(rref(B + [b])) > len(B): B.append(b)
    assert len(B) == m + 1
    nb = []
    # dual vectors: v_i . B_j = delta_ij (j = 1..m), v_i . ONE = 0   (solve over F_2 by brute elimination)
    eqs = [B[0]] + B[1:]
    duals = []
    for i in range(1, m + 1):
        rhs = [0] + [1 if j == i else 0 for j in range(1, m + 1)]
        # solve A v = rhs with A rows = eqs (as bitmasks over n unknowns)
        rows_ = [(e, r) for e, r in zip(eqs, rhs)]
        piv = []
        for k in range(len(rows_)):
            e, r = rows_[k]
            for (pe, pr, pb) in piv:
                if e >> pb & 1: e ^= pe; r ^= pr
            assert e
            pb = e.bit_length() - 1
            piv = [((qe ^ e, qr ^ r, qb) if (qe >> pb & 1) else (qe, qr, qb)) for (qe, qr, qb) in piv]
            piv.append((e, r, pb))
        v = 0
        for (pe, pr, pb) in piv:
            if pr: v |= 1 << pb
        for e, r in zip(eqs, rhs): assert dot(v, e) == r
        duals.append(v)
    for mask in range(1, 1 << m):
        v = 0
        for i in range(m):
            if mask >> i & 1: v ^= duals[i]
        C0 = [B[0]] + [B[i] for i in range(1, m + 1) if not (mask >> (i - 1) & 1)]
        ones = [B[i] for i in range(1, m + 1) if (mask >> (i - 1) & 1)]
        b = ones[0]
        C0 += [ones[0] ^ o for o in ones[1:]]
        for D in (rref(C0 + [v]), rref(C0 + [v ^ b])):
            if D not in idx: idx[D] = len(codes); codes.append(D)
            nb.append(idx[D])
    nbrs.append(nb)
N = len(codes)
expect = 1
for i in range(1, m + 1): expect *= (2 ** i + 1)
print(f"m={m}: frames {N} (expected {expect}), degree {len(nbrs[0])} (expected {2*(2**m-1)})", file=sys.stderr)
assert N == expect
for i in range(N):
    for j in nbrs[i]: assert i in nbrs[j]
def is_pair(C):
    w2 = [w for w in span(C) if bin(w).count('1') == 2]
    cov = 0
    for w in w2:
        if cov & w: return None
        cov |= w
    return w2 if cov == ONE else None
out = open(sys.argv[2], 'w')
out.write(f"{N} {m} {d} {2*(2**m-1)}\n")
npair = 0
for i, C in enumerate(codes):
    pr = is_pair(C)
    L = rref([w for w in span(C) if not (w >> INF & 1)])
    assert len(L) == m
    rows = []
    if pr is not None:
        npair += 1
        for w in pr:
            if w >> INF & 1: continue
            a, b = [k for k in range(n) if w >> k & 1]
            r = [0] * d; r[a] = 1; r[b] = -1; rows.append(r)
    else:
        w2 = [w for w in span(C) if bin(w).count('1') == 2]   # pair components (weight-2 words)
        cov = 0
        for w in w2: cov |= w
        for w in w2:
            if w >> INF & 1: continue
            a, b = [k for k in range(n) if w >> k & 1]
            r = [0] * d; r[a] = 1; r[b] = -1; rows.append(r)
        restL = rref([w for w in span(C) if not (w >> INF & 1) and not (w & cov)])
        for w in restL:
            pts = [k for k in range(d) if w >> k & 1]
            r = [0] * d
            for t, p in enumerate(pts): r[p] = 1 if t % 2 == 0 else -1   # signed: row sums to 0, like e_a - e_b
            rows.append(r)
    assert len(rows) == m
    # LIFT CHECK: reductions mod 2 of the rows are m independent words of L (so they form a basis of L)
    red = [sum(((x & 1) << k) for k, x in enumerate(r)) for r in rows]
    Lw = set(span(L))
    assert all(w in Lw for w in red) and len(rref(red)) == m
    out.write(("1 " if pr is not None else "0 ") + " ".join(str(x) for r in rows for x in r) + " " +
              " ".join(map(str, nbrs[i])) + "\n")
out.close()
print(f"pair frames {npair}", file=sys.stderr)
