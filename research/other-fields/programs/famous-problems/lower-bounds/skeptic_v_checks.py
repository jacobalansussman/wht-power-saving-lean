"""v-iso-lower independent checks (small; < 50 MB, seconds).

(A) NEW: Bluestein analogue for the Walsh-Hadamard transform.
    Over the integers, |x XOR y| = |x| + |y| - 2 |x AND y|, so 2 (x.y) = |x| + |y| - |x XOR y| and
    (-1)^(x.y) = i^|x| * i^|y| * i^(-|x XOR y|).
    Hence H v = D * conv_xor(D v, c) with D = diag(i^|y|), c(z) = i^(-|z|):
    the WHT is two diagonal scalings plus ONE dyadic (XOR) convolution with a fixed vector.
    With the dyadic convolution theorem conv_xor(a, b) = H(Ha o Hb)/N this gives
    WHT  <=  2N + 3 * XORCONV(N)   and   XORCONV(N) <= 3 WHT(N) + 2N,
    so the WHT exponent equals the dyadic-convolution exponent (EQ), the WHT twin of
    DFT <-> polynomial multiplication.
(B) Independent re-check of Bluestein (DFT via one linear convolution with a fixed chirp) and of the
    convolution theorem, all lengths 1..40.
(C) Linear independence of f -> P^(tensor f) for distinct nonzero matrices: generic x, y give
    sum_i c_i (y^T P_i x)^f = 0 for all f, a Vandermonde system with distinct nonzero nodes.
(D) The radix-2 butterfly (FFT graph for H) is not a 2-superconcentrator, and the matching 2x2
    minor of H vanishes, for n = 3, 4, 5.
(E) Morgenstern potential, both gate models: addition gate factor <= 2, scalar gate factor <= |c|;
    a fused gate a u + b v has factor |a| + |b| (<= 2c), so the base is max(2,c) in OpenAI's model
    (u+v, u-v, lambda u each one gate) and 2c in the fused model (Raz's log_{2c}).
"""
import itertools
import numpy as np

rng = np.random.default_rng(1)


def wht_matrix(n):
    H = np.array([[1.0]])
    h = np.array([[1.0, 1.0], [1.0, -1.0]])
    for _ in range(n):
        H = np.kron(H, h)
    return H


def popcount(a):
    return bin(a).count("1")


def xor_conv(a, b):
    N = len(a)
    out = np.zeros(N, dtype=complex)
    for x in range(N):
        for y in range(N):
            out[x ^ y] += a[x] * b[y]
    return out


print("(A) WHT = diag * XOR-convolution with fixed c(z)=i^-|z| * diag")
for n in range(1, 8):
    N = 2 ** n
    H = wht_matrix(n)
    v = rng.standard_normal(N) + 1j * rng.standard_normal(N)
    D = np.array([1j ** popcount(y) for y in range(N)])
    c = np.array([(1j) ** (-popcount(z)) for z in range(N)])
    # (Hv)_y = i^|y| sum_x (i^|x| v_x) c(x XOR y)
    w = D * v
    conv = xor_conv(w, c)  # conv[y] = sum_x w_x c(x ^ y) since (x ^ (x^y)) = y
    lhs = H @ v
    rhs = D * conv
    # also the convolution theorem for XOR convolution
    a = rng.standard_normal(N); b = rng.standard_normal(N)
    ct = H @ ((H @ a) * (H @ b)) / N
    print(f"   n={n} N={N}: max|Hv - D*xorconv(Dv,c)| = {np.max(np.abs(lhs - rhs)):.2e};"
          f" max|xorconv - H(Ha.Hb)/N| = {np.max(np.abs(xor_conv(a, b) - ct)):.2e}")

print("(B) Bluestein and convolution theorem, n = 1..40")
worst_b = worst_c = 0.0
for n in range(1, 41):
    F = np.exp(2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n)
    x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    j = np.arange(n)
    # jk = (j^2 + k^2 - (k-j)^2)/2 ; zeta^{jk} = w(j) w(k) / w(k-j) with w(t) = zeta^{t^2/2}
    w = lambda t: np.exp(1j * np.pi * (t.astype(float) ** 2) / n)
    a = x * w(j)                                  # n scalings
    m = np.arange(-(n - 1), n)                    # chirp indices -(n-1)..(n-1)
    chirp = 1.0 / w(m)                            # fixed polynomial of degree 2n-2
    prod = np.convolve(a, chirp)                  # one product: variable deg n-1 times fixed deg 2n-2
    X = w(j) * prod[n - 1: 2 * n - 1]             # middle coefficients, n scalings
    worst_b = max(worst_b, np.max(np.abs(X - F @ x)) / max(1, np.max(np.abs(F @ x))))
    # convolution theorem with length L >= 2n-1 (power of two <= 4n)
    p = rng.standard_normal(n); q = rng.standard_normal(n)
    L = 1
    while L < 2 * n - 1:
        L *= 2
    FL = np.exp(-2j * np.pi * np.outer(np.arange(L), np.arange(L)) / L)
    pp = np.zeros(L); pp[:n] = p
    qq = np.zeros(L); qq[:n] = q
    cc = np.conj(FL) @ ((FL @ pp) * (FL @ qq)) / L
    worst_c = max(worst_c, np.max(np.abs(cc[:2 * n - 1] - np.convolve(p, q))))
    assert L <= 4 * n
print(f"   Bluestein worst relative error {worst_b:.2e}; convolution theorem worst error {worst_c:.2e}; L <= 4n always")

print("(C) linear independence of f -> P^(x)f via Vandermonde at generic x, y")
for q, k in [(2, 12), (3, 9)]:
    Ps = []
    while len(Ps) < k:
        P = rng.integers(-1, 2, size=(q, q)).astype(float)
        if np.any(P) and all(not np.array_equal(P, Q) for Q in Ps):
            Ps.append(P)
    x = rng.standard_normal(q); y = rng.standard_normal(q)
    lam = np.array([y @ P @ x for P in Ps])
    V = np.array([[l ** f for l in lam] for f in range(1, k + 1)])
    print(f"   q={q}, {k} distinct nonzero matrices: nodes distinct & nonzero: "
          f"{len(set(np.round(lam, 12))) == k and np.all(np.abs(lam) > 1e-9)}; Vandermonde rank {np.linalg.matrix_rank(V)} / {k}")

print("(D) butterfly graph of the FFT for H_n is not a 2-superconcentrator")
for n in range(3, 6):
    N = 2 ** n
    H = wht_matrix(n)
    # path from input x to output y in the decimation network: after stage j the position holds the
    # top j bits (processed first) from y and the remaining bits from x.  Bits processed in order 0..n-1.
    def node(x, y, j):
        mask = (1 << j) - 1
        return (j, (y & mask) | (x & ~mask))
    x1 = 0b0; x2 = 0b1                      # differ in bit 0 (processed first)
    y1 = 0b0; y2 = 1 << (n - 1)             # differ in bit n-1 (processed last)
    collide = []
    for (a, b), (c_, d) in [((x1, y1), (x2, y2)), ((x1, y2), (x2, y1))]:
        p1 = {node(a, b, j) for j in range(1, n)}
        p2 = {node(c_, d, j) for j in range(1, n)}
        collide.append(len(p1 & p2) > 0)
    minor = H[y1, x1] * H[y2, x2] - H[y1, x2] * H[y2, x1]
    print(f"   n={n}: both pairings share an internal vertex: {all(collide)}; 2x2 minor of H on those rows/cols = {minor}")

print("(E) Morgenstern potential factors")
for trial in range(3):
    n = 4
    rows = np.eye(n)
    best = 1.0
    worst_add = worst_scal = 0.0
    for g in range(30):
        i, j = rng.integers(0, len(rows), 2)
        before = max(abs(np.linalg.det(rows[list(S)])) for S in itertools.combinations(range(len(rows)), n))
        if rng.random() < 0.5:
            new = rows[i] + rows[j]
            kind = "add"
        else:
            lam = rng.uniform(-3, 3)
            new = lam * rows[i]
            kind = "scal"
        rows = np.vstack([rows, new])
        after = max(abs(np.linalg.det(rows[list(S)])) for S in itertools.combinations(range(len(rows)), n)) if len(rows) <= 12 else None
        if after is None:
            break
        ratio = after / before
        if kind == "add":
            worst_add = max(worst_add, ratio)
        else:
            worst_scal = max(worst_scal, ratio / max(1.0, abs(lam)))
    print(f"   trial {trial}: worst addition-gate factor {worst_add:.3f} (<= 2); worst scalar-gate factor / max(1,|lam|) {worst_scal:.3f} (<= 1)")
