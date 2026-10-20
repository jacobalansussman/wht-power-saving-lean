"""toy_d.py (u-ours, wht43/copywall) -- can a LINEAR total (a functional of the root vector: the 'natural family' of
wht35/e8/s-copies section 4, the family of the hub's coordinate totals) on P112 = E8 minus one octad stand at a frame
of dimension 7 or 6?  For each subspace F of F_2^9 of dimension d, the P-ports outside F must all be zeros of the
functional, so their roots must have rank < 8.  Rank via coordinates on the removed octad (an orthogonal basis of
R^8): small exact integers.  usage: python3 -B toy_d.py [octad index, default 0]"""
import os
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
import sys, itertools, collections
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "../e8-labels"))
import e8
import numpy as np
from fractions import Fraction as Fr
V = e8.V; PORTS = e8.PORTS
GRAM = np.array(e8.GRAM, dtype=np.int64); ORTH = (GRAM == 0)
def par(x): return bin(x).count("1") & 1
OM = [sum(1 << T for T in range(V) if ORTH[S, T]) for S in range(V)]
def cliques(adj, size, uni):
    out = []
    def rec(cur, cand):
        if len(cur) == size: out.append(tuple(cur)); return
        c = cand
        while c:
            low = c & -c; i = low.bit_length() - 1; c ^= low
            if bin(cand >> i).count("1") + len(cur) < size: break
            rec(cur + [i], cand & adj[i] & ~((low << 1) - 1))
    rec([], uni); return out
OCT = cliques(OM, 8, (1 << V) - 1)
O = OCT[int(sys.argv[1]) if len(sys.argv) > 1 else 0]
P = [T for T in range(V) if T not in O]
Pmask = np.zeros(V, bool); Pmask[P] = True
X = GRAM[:, list(O)]                      # 120 x 8 coordinates (times 2/|w|^2) on the octad basis
assert np.linalg.matrix_rank(X.astype(float)) == 8
def frank(rows):
    """exact rank over Q by fraction-free elimination on small ints"""
    A = [list(map(int, r)) for r in rows]; rk = 0; ncol = 8
    for c in range(ncol):
        piv = next((i for i in range(rk, len(A)) if A[i][c]), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][c]:
                a, b = A[rk][c], A[i][c]
                A[i] = [a * y - b * x for x, y in zip(A[rk], A[i])]
        rk += 1
        if rk == 8: break
    return rk
P0 = np.array([[par(PORTS[T] & c) == 0 and Pmask[T] for T in range(V)] for c in range(512)], dtype=bool)
P0pk = np.packbits(P0, axis=1)
def subspaces(k, n=9):
    for piv in itertools.combinations(range(n - 1, -1, -1), k):
        free = [[j for j in range(p) if j not in piv] for p in piv]
        f = sum(len(x) for x in free); idx = np.arange(1 << f, dtype=np.int64); cols = []; sh = 0
        for p, fr in zip(piv, free):
            col = np.full(1 << f, 1 << p, dtype=np.int64)
            for j in fr: col |= ((idx >> sh) & 1) << j; sh += 1
            cols.append(col)
        yield np.stack(cols, axis=1)
for d in (8, 7, 6):
    k = 9 - d; seen = set()
    for bas in subspaces(k):
        m = P0pk[bas[:, 0]]
        for i in range(1, k): m = m & P0pk[bas[:, i]]
        for row in np.unique(m, axis=0): seen.add(row.tobytes())
    cnt = collections.Counter(); low = []
    for key in seen:
        inside = np.unpackbits(np.frombuffer(key, np.uint8))[:V].astype(bool)
        outside = [T for T in P if not inside[T]]
        fl = np.linalg.matrix_rank(X[outside].astype(float)) if outside else 0
        if fl < 8:
            fl = frank(X[outside]); low.append((fl, int(inside.sum())))
        cnt[fl] += 1
    print("d = %d: %d distinct P-port sets; rank of the roots outside the frame: %s; frames admitting a linear total "
          "(rank < 8): %d, their (rank, P-ports inside): %s" % (d, len(seen), sorted(cnt.items()), len(low),
          sorted(collections.Counter(low).items())), flush=True)
