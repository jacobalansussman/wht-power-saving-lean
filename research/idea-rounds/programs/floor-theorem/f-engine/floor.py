#!/usr/bin/python3
"""floor.py (agent f-engine, 10 Oct 2026).  Checks for the floor theorem "Delta certificate".
Matrices typed from the Lean definitions:
  dir z = ca*1 + cb*shift z, ca=(1+i)/2, cb=(1-i)/2      OAI/.../BinaryFrames.lean:67-69,152
  shift z i j = 1 iff i+z=j                               BinaryFrames.lean:137-138
  permMat G x y = 1 iff G.pi x = y                        Work/GFrame/Engine/Letters.lean:39
  phaseMat z c = diag(tint(z.x+c)), tint 0=1, tint 1=i    Letters.lean:42
  block = dir z1 ... dir zk, z independent, 1<=k<m        Work/Block/Split.lean:125, Word.lean:30,55-57
All entries are Gaussian numbers with a power of 2 below, so complex128 arithmetic is exact here;
ranks are computed exactly over Q(i) with Fractions (part A) or modulo the prime 65537 with i -> 256
(part B; a rank modulo a prime is a lower bound for the true rank) and cross-checked by SVD.
NOT Lean.  Nobody has read this.
"""
import sys, itertools, random, math
from fractions import Fraction as Fr
import numpy as np
np.seterr(all="ignore")

def dot(x, y): return bin(x & y).count("1") & 1
def shift(m, z):
    n = 1 << m; a = np.zeros((n, n), complex)
    for i in range(n): a[i, i ^ z] = 1
    return a
def dirm(m, z): return (1 + 1j) / 2 * np.eye(1 << m) + (1 - 1j) / 2 * shift(m, z)
def phase(m, z, c): return np.diag([1j if (dot(z, x) + c) & 1 else 1 for x in range(1 << m)]).astype(complex)
def lin(m, cols, x):
    y = 0
    for i in range(m):
        if x >> i & 1: y ^= cols[i]
    return y
def perm(m, cols):
    n = 1 << m; a = np.zeros((n, n), complex)
    for x in range(n): a[x, lin(m, cols, x)] = 1
    return a
def indep(m, zs):
    span = {0}
    for z in zs:
        if z in span: return False
        span |= {s ^ z for s in span}
    return True
def block(m, zs):
    M = np.eye(1 << m, dtype=complex)
    for z in zs: M = M @ dirm(m, z)
    return M
def kron_f(M, f):
    out = np.ones((1, 1), complex)
    for _ in range(f): out = np.kron(out, M)
    return out

def rank_exact(M):
    """exact rank over Q(i) of a complex array whose entries are exact dyadic Gaussian numbers"""
    rows = [[(Fr(v.real), Fr(v.imag)) for v in r] for r in M]
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        p = None
        for r in range(rk, len(rows)):
            if rows[r][c] != (0, 0): p = r; break
        if p is None: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        a, b = rows[rk][c]; d = a * a + b * b; ia, ib = a / d, -b / d
        for r in range(rk + 1, len(rows)):
            u, v = rows[r][c]
            if (u, v) == (0, 0): continue
            fa, fb = u * ia - v * ib, u * ib + v * ia
            rows[r] = [(x - (fa * s - fb * t), y - (fa * t + fb * s)) for (x, y), (s, t) in zip(rows[r], rows[rk])]
        rk += 1
        if rk == len(rows): break
    return rk

P = 65537; IM = 256           # 256^2 = -1 mod 65537
def rank_modp(M):
    A = np.round(M * 2.0 ** 30)
    assert np.array_equal(A, M * 2.0 ** 30), "entries are not dyadic at 2^-30"
    Z = (A.real.astype(np.int64) % P + IM * (A.imag.astype(np.int64) % P)) % P
    Z = Z.copy(); rk = 0; nr, nc = Z.shape
    for c in range(nc):
        nz = np.nonzero(Z[rk:, c])[0]
        if len(nz) == 0: continue
        p = rk + nz[0]
        if p != rk: Z[[rk, p]] = Z[[p, rk]]
        inv = pow(int(Z[rk, c]), P - 2, P)
        Z[rk] = Z[rk] * inv % P
        for r in range(rk + 1, nr):
            if Z[r, c]: Z[r] = (Z[r] - Z[r, c] * Z[rk]) % P
        rk += 1
        if rk == nr: break
    return rk

def diag_basis(m, j=1):
    """columns e_X (x) e_X, X an address tuple of j columns: the subspace Delta^(j) of V^(2j)"""
    n = (1 << m) ** j; B = np.zeros((n * n, n), complex)
    for x in range(n): B[x * n + x, x] = 1
    return B

def partA(m, sample=None):
    n = 1 << m; D = diag_basis(m); ok = True
    # (1) every free letter keeps Delta
    frees = [shift(m, z) for z in range(1, n)] + [phase(m, z, c) for z in range(n) for c in (0, 1)]
    gl = [cols for cols in itertools.product(range(1, n), repeat=m) if indep(m, cols)]
    if sample: random.seed(1); gl = random.sample(gl, sample)
    frees += [perm(m, cols) for cols in gl]
    offd = [x * n + y for x in range(n) for y in range(n) if x != y]
    for b in frees:
        img = np.kron(b, b) @ D
        if np.any(img[offd] != 0) or abs(np.linalg.det(b)) == 0: ok = False
    print(f"m={m}: {len(frees)} free letters ({n-1} shifts, {2*n} phases, {len(gl)} address maps): "
          f"each maps Delta into Delta: {ok}")
    # (2) weight of every block
    for k in range(1, m + 1):
        tuples = [zs for zs in itertools.product(range(1, n), repeat=k) if indep(m, zs)]
        if sample and len(tuples) > sample: random.seed(k); tuples = random.sample(tuples, sample)
        vals = set()
        for zs in tuples:
            b = block(m, zs)
            r = rank_exact(np.hstack([D, np.kron(b, b) @ D]).T)
            vals.add(r - n)
        name = "kernel-rank" if k == m else "block"
        print(f"   rank {k} ({name}; {len(tuples)} tuples of independent directions): "
              f"dim(Delta + (g(x)g)Delta) - dim Delta = {sorted(vals)}   predicted 2^m - 2^(m-k) = {n - (n >> k)}")
    K = np.ones((1, 1), complex)
    for _ in range(m): K = np.kron(K, np.array([[(1 + 1j) / 2, (1 - 1j) / 2], [(1 - 1j) / 2, (1 + 1j) / 2]]))
    print("   kernel (digitProd littleC) == dir e1 ... dir em :", np.array_equal(K, block(m, [1 << i for i in range(m)])))

def partA_higher(m):
    """columns f = 4 (j = 2): weight 2^(2m) - (2^(m-k))^2"""
    n = 1 << m; D = diag_basis(m, 2)
    for k in range(1, m + 1):
        b = block(m, [1 << i for i in range(k)]); b2 = np.kron(b, b)
        r = rank_modp(np.hstack([D, np.kron(b2, b2) @ D]).T)
        r2 = np.linalg.matrix_rank(np.hstack([D, np.kron(b2, b2) @ D]))
        print(f"m={m}, 4 columns, rank {k}: weight = {r - n*n} (mod p), {r2 - n*n} (SVD); predicted {n*n - (n >> k)**2}")

# ---------- part B: the potential along words, m = 3, two columns ----------
def run_word(m, nrole, live, word, f=2, track=True):
    n = (1 << m) ** f; W = len(live)
    T = np.zeros((nrole, n, W * n), complex)
    for l, r in enumerate(live): T[r, :, l * n:(l + 1) * n] = np.eye(n)
    nn = 1 << m
    ycols = [l * n + x * nn + x for l in range(W) for x in range(nn)] if f == 2 else None
    offd = [x * nn + y for x in range(nn) for y in range(nn) if x != y] if f == 2 else None
    def phi():
        M = np.vstack([T[r][np.ix_(offd, ycols)] for r in range(nrole)])
        a, b = rank_modp(M), int(np.linalg.matrix_rank(M))
        assert a == b, (a, b)
        return a
    hist = [phi()] if (track and f == 2) else []
    viol = 0; paid = 0
    for let in word:
        if let[0] == "add":
            T = np.einsum("rj,jab->rab", np.array(let[1], dtype=complex), T); cap = 0
        elif let[0] == "free":
            T[let[1]] = kron_f(let[2], f) @ T[let[1]]; cap = 0
        else:
            _, r, zs = let
            T[r] = kron_f(block(m, zs), f) @ T[r]; cap = nn - (nn >> len(zs)); paid += cap
        if track and f == 2:
            hist.append(phi())
            if hist[-1] > hist[-2] + cap: viol += 1
    K = kron_f(block(m, [1 << i for i in range(m)]), f)
    good = all(np.array_equal(T[r][:, l * n:(l + 1) * n], K) and
               all(not np.any(T[r][:, q * n:(q + 1) * n]) for q in range(W) if q != l)
               for l, r in enumerate(live))
    return good, hist, viol, paid

def partB():
    m = 3; e1, e2, e3 = 1, 2, 4
    I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    plain = [("blk", 0, (e1, e2)), ("blk", 0, (e3,)), ("blk", 1, (e1, e2)), ("blk", 1, (e3,))]
    # "sum trick": live arrays 0 and 1, scratch 2.  s = a0 + a1 ; kernel on s and on a0 ; a1 := s - a0.
    # Array 1 is never the target of a block.
    trick = [("add", [[1, 0, 0], [0, 1, 0], [1, 1, 0]]),
             ("blk", 2, (e1, e2)), ("blk", 2, (e3,)), ("blk", 0, (e1, e2)), ("blk", 0, (e3,)),
             ("add", [[1, 0, 0], [-1, 0, 1], [0, 0, 1]])]
    # a wrong word for contrast: one block of rank 2 per array and free letters only
    wrong = [("blk", 0, (e1, e2)), ("free", 0, phase(m, 4, 0)), ("free", 0, shift(m, 4)),
             ("blk", 1, (e1, e2)), ("free", 1, phase(m, 4, 0)), ("free", 1, shift(m, 4))]
    for name, w in (("plain (ranks 2,1 per array)", plain), ("sum trick (array 1 gets no block)", trick),
                    ("one rank-2 block per array + free letters", wrong)):
        res = []
        for f in (0, 1, 2, 3):
            good, hist, viol, paid = run_word(m, 3, [0, 1], w, f=f, track=(f == 2))
            res.append(f"f={f}:{'right' if good else 'WRONG'}")
            if f == 2: keep = (hist, viol, paid)
        hist, viol, paid = keep
        print(f"{name}: {' '.join(res)} | potential along the word {hist} | sum of block weights {paid} | "
              f"needed 2*(8-1)=14 | steps breaking the lemma: {viol}")
    # random words: the lemma's two inequalities are claimed for EVERY word, right or wrong
    random.seed(20261010)
    n = 8; gl = [cols for cols in itertools.product(range(1, n), repeat=3) if indep(3, cols)]
    coef = [0, 0, 0, 1, -1, Fr(1, 2), Fr(-1, 2), 2, Fr(3, 4)]
    tot = 0; viols = 0; rises = 0; maxphi = 0
    for trial in range(120):
        word = []
        for _ in range(24):
            t = random.random()
            if t < 0.35:
                word.append(("add", [[float(random.choice(coef)) + (1.0 if i == j else 0.0) * (random.random() < 0.7)
                                      for j in range(3)] for i in range(3)]))
            elif t < 0.6:
                kind = random.randrange(3); r = random.randrange(3)
                Mx = shift(3, random.randrange(1, 8)) if kind == 0 else \
                    (phase(3, random.randrange(8), random.randrange(2)) if kind == 1 else perm(3, random.choice(gl)))
                word.append(("free", r, Mx))
            else:
                k = random.choice((1, 2))
                while True:
                    zs = tuple(random.randrange(1, 8) for _ in range(k))
                    if indep(3, zs): break
                word.append(("blk", random.randrange(3), zs))
        good, hist, viol, paid = run_word(3, 3, [0, 1], word, f=2)
        tot += len(word); viols += viol; maxphi = max(maxphi, max(hist))
        rises += sum(1 for a, b in zip(hist, hist[1:]) if b > a)
    print(f"random words: 120 words, {tot} letters, potential rose at {rises} letters, largest potential seen {maxphi} "
          f"(ceiling 3*56=168), letters breaking the lemma: {viols}")

def partC():
    print("m | blocks per live array at least (2^m-1)/(2^m-2) | z_0(m) = min_k ln(w(m)/w(k))/ln(m/k), w(k)=1-2^-k | attained at k | moves per array at least 2(1-2^-m)")
    for m in (3, 4, 5, 6, 8, 10, 12, 16, 24, 45, 110):
        best = None
        for k in range(1, m):
            num = math.log1p(-2.0 ** -m) - math.log1p(-2.0 ** -k)
            z = num / math.log(m / k)
            if best is None or z < best[0]: best = (z, k)
        print(f"{m:4d} | 1 + {1.0/(2.0**m-2):.6e} | {best[0]:.6e} | k={best[1]} | {2*(1-2.0**-m):.6f}")
    # what the table of wht35/wider/w-bounds.md (block size 3; argued there, not by me) would give
    z1 = math.log(3 / (41 / 35)) / math.log(3); z2 = math.log(3 / (1307 / 630)) / math.log(3 / 2)
    print(f"m=3 with w-bounds' table n1*41/35 + n2*1307/630 >= 3W: z >= min({z1:.4f}, {z2:.4f}) = {min(z1,z2):.4f}; "
          f"blocks per array >= {3*630/1307:.4f}")
    # the recurrence at exactly one block of rank m-1 per array, and at fewer blocks than arrays
    for m, lam in ((3, 1.0), (3, 0.9), (110, 1.0)):
        def tt(k):
            chain = []
            while k >= 3 * m: chain.append(k); k = (m - 1) * (k // m)
            v = 1.0
            for _ in chain: v = 1.0 + lam * v
            return v
        ks = [10 ** 3, 10 ** 6, 10 ** 12, 10 ** 24]
        print(f"t(k) = 1 + {lam}*t((m-1)*floor(k/m)), m={m}: " + ", ".join(f"t(1e{len(str(k))-1})={tt(k):.2f}" for k in ks)
              + (f" ; ln k/ln(m/(m-1)) at 1e24 = {math.log(1e24)/math.log(m/(m-1)):.1f}" if lam == 1.0 else f" ; limit 1/(1-lam) = {1/(1-lam):.1f}"))

if __name__ == "__main__":
    part = sys.argv[1]
    if part == "A": partA(2); partA(3); partA(4, sample=40); partA_higher(2)
    if part == "B": partB()
    if part == "C": partC()
