"""vland.py (v-fp-build): Lemma 3.2 (landing) with my own constructive chart (vfg.py arithmetic only).
For random walker g, random stretch w (product of transvections with nonzero x-part and random free letters), and a
random target X with d(g^-1 L0, X) = level(w): build mu in P with pos(w mu g) = X by my own normal form
(Levi straightening of the x-projection, then a symmetric phase), and check. Also check the converse: when
d(g^-1 L0, X) != level(w), no mu in P works (distance invariant), on the same random data."""
import sys, random
from vfg import FG
def randP(S, rng):
    m = S.m
    while True:  # Levi: random invertible A on x, A^-T on z
        cols = [rng.randrange(1, 1 << m) for _ in range(m)]
        if S.rank(cols) == m: break
    A = tuple(cols) + tuple(0 for _ in range(m))
    # build g: (a|b) -> (A a | A^-T b + C a); first the Levi part
    Am = [cols[j] for j in range(m)]            # column j of A
    # A^-T: compute inverse of A over F2
    def matinv(colsA):
        N = m; rows = [sum(((colsA[j] >> i) & 1) << j for j in range(N)) | (1 << (N + i)) for i in range(N)]
        for c in range(N):
            p = next(r for r in range(c, N) if (rows[r] >> c) & 1); rows[c], rows[p] = rows[p], rows[c]
            for r in range(N):
                if r != c and (rows[r] >> c) & 1: rows[r] ^= rows[c]
        inv_rows = [rows[i] >> N for i in range(N)]
        return [sum(((inv_rows[i] >> j) & 1) << i for i in range(N)) for j in range(N)]  # columns of A^-1
    Ainv = matinv(Am)
    # A^-T column j = row j of A^-1
    AinvT = [sum(((Ainv[i] >> j) & 1) << i for i in range(m)) for j in range(m)]
    levi = tuple(Am[j] for j in range(m)) + tuple(AinvT[j] << m for j in range(m))
    # random symmetric phase: (a|b) -> (a | b + Sym a)
    Sym = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i, m):
            Sym[i][j] = Sym[j][i] = rng.randrange(2)
    ph = tuple((1 << j) | (sum(Sym[i][j] << i for i in range(m)) << m) for j in range(m)) + tuple(1 << (m + j) for j in range(m))
    return S.mul(ph, levi)
def to_std(S, Y):
    """my normal form: p in P with p Y = L_d = span{(e_1|0)..(e_d|0), (0|e_{d+1})..(0|e_m)}"""
    m = S.m; M = (1 << m) - 1
    xs = S.rref([v & M for v in Y]); d = len(xs)
    # Levi A with A(xs basis) = e_1..e_d : extend xs to a basis, A = inverse of [xs | ext]
    basis = list(xs)
    for e in range(m):
        if S.rank(basis + [1 << e]) > len(basis): basis.append(1 << e)
    # B has columns basis; A = B^-1
    N = m; rows = [sum(((basis[j] >> i) & 1) << j for j in range(N)) | (1 << (N + i)) for i in range(N)]
    for c in range(N):
        p = next(r for r in range(c, N) if (rows[r] >> c) & 1); rows[c], rows[p] = rows[p], rows[c]
        for r in range(N):
            if r != c and (rows[r] >> c) & 1: rows[r] ^= rows[c]
    inv_rows = [rows[i] >> N for i in range(N)]
    A = [sum(((inv_rows[i] >> j) & 1) << i for i in range(N)) for j in range(N)]   # columns of B^-1
    AT_rows = basis  # A^-T = B^T: column j of B^T = row j of B
    BT = [sum(((basis[i] >> j) & 1) << i for i in range(m)) for j in range(m)]
    levi = tuple(A) + tuple(BT[j] << m for j in range(m))
    Y1 = S.img(levi, Y)
    # now Y1 = span{(e_i | s_i), i < d} + span{(0|e_j), j >= d}; read s_i restricted to first d coords
    span = [0]
    for v in Y1: span = span + [x ^ v for x in span]
    srow = {}
    for v in span:
        a = v & M
        if a and a & (a - 1) == 0: srow[a.bit_length() - 1] = v
    Sym = [[0] * m for _ in range(m)]
    for i in range(d):
        v = srow[i]
        assert v & M == 1 << i
        s = (v >> m) & ((1 << d) - 1)
        for j in range(d): Sym[j][i] = (s >> j) & 1
    ph = tuple((1 << j) | (sum(Sym[i][j] << i for i in range(m)) << m) for j in range(m)) + tuple(1 << (m + j) for j in range(m))
    p = S.mul(ph, levi)
    std = S.rref([1 << i for i in range(d)] + [1 << (m + j) for j in range(d, m)])
    assert S.inP(p) and S.img(p, Y) == std, "normal form failed"
    return p
rng = random.Random(77)
for m in (3, 4, 5, 6):
    S = FG(m); ok = 0; neg_ok = 0; neg = 0
    xvecs = [u for u in range(1, 1 << (2 * m)) if u & ((1 << m) - 1)]
    for t in range(200):
        g = S.I
        for _ in range(rng.randrange(1, 2 * m + 2)): g = S.mul(S.tv(rng.choice(xvecs)), S.mul(randP(S, rng), g))
        w = S.I
        for _ in range(rng.randrange(1, m + 2)): w = S.mul(S.tv(rng.choice(xvecs)), S.mul(randP(S, rng), w))
        lw = S.level(w)
        # random X at distance lw from pos(g)
        pg = S.pos(g)
        Ld = S.rref([1 << i for i in range(lw)] + [1 << (m + j) for j in range(lw, m)])
        Y = S.img(randP(S, rng), Ld)           # random point at distance lw from L0
        X = S.img(S.inv(g), Y)                 # at distance lw from pos(g) = g^-1 L0
        assert S.dist(pg, X) == lw
        # need mu with mu(gX) = w^-1 L0
        gX = S.img(g, X); wL = S.img(S.inv(w), S.L0)
        mu = S.mul(S.inv(to_std(S, wL)), to_std(S, gX))
        if S.inP(mu) and S.pos(S.mul(w, S.mul(mu, g))) == X: ok += 1
        # negative: X' at another distance: then pos(w mu g) != X' for random mu (invariant)
        h = S.I
        for _ in range(3): h = S.mul(S.tv(rng.choice(xvecs)), h)
        X2 = S.pos(h)
        if S.dist(pg, X2) != lw:
            neg += 1
            mu2 = randP(S, rng)
            if S.dist(pg, S.pos(S.mul(w, S.mul(mu2, g)))) == lw: neg_ok += 1
    print(f"m={m}: landing chart correct in {ok}/200; invariant d(pos g, pos(w mu g)) = level(w) for random mu: {neg_ok}/{neg}")
