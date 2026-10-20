"""fastP.py (w3-b5): exact, uniform sampling of {g in P : g X = Y} (P = Stab(L0)) without rejection.
s_X in P maps X to the canonical C_k = E_k (+) Z(E_k^perp) (k = d(L0, X)); the coset is s_Y^-1 Stab_P(C_k) s_X, and
Stab_P(C_k) = {unip(S) levi(A): A E_k = E_k, S E_k in E_k^perp} is sampled uniformly."""
import random
from sp2 import M, MASKX, L0, levi, unip, mul, inv, image, rank, span_all, rref, dist, rand_gl


def inv5(cols):
    """inverse of an invertible m x m F2 matrix given by column images (list of m ints)"""
    n = M
    rows = []
    for r in range(n):
        row = 0
        for i in range(n):
            if (cols[i] >> r) & 1:
                row |= 1 << i
        rows.append([row, 1 << r])
    piv = 0
    for c in range(n):
        pr = next(r for r in range(piv, n) if (rows[r][0] >> c) & 1)
        rows[piv], rows[pr] = rows[pr], rows[piv]
        for r in range(n):
            if r != piv and (rows[r][0] >> c) & 1:
                rows[r][0] ^= rows[piv][0]; rows[r][1] ^= rows[piv][1]
        piv += 1
    out = []
    for i in range(n):
        col = 0
        for r in range(n):
            if (rows[r][1] >> i) & 1:
                col |= 1 << r
        out.append(col)
    return out


def canon(k):
    return rref([1 << i for i in range(k)] + [1 << (M + i) for i in range(k, M)])


def s_of(X):
    """s in P with s X = canon(k)"""
    S_ = span_all(X)
    U = []
    for v in sorted(S_):
        a = v & MASKX
        if a and rank(U + [a]) > len(U):
            U.append(a)
    k = len(U)
    cols = list(U)
    for i in range(M):
        if rank(cols + [1 << i]) > len(cols):
            cols.append(1 << i)
    A = inv5(cols)                      # A u_i = e_i
    lv = levi(A)
    Xp = image(lv, X)
    Sp_ = span_all(Xp)
    S = [0] * M
    for i in range(k):
        v = next(v for v in Sp_ if (v & MASKX) == (1 << i))
        b = v >> M
        S[i] = b & ((1 << k) - 1)
    g = mul(unip(S), lv)
    assert image(g, X) == canon(k), "canon failed"
    return g, k


def rand_stab(rng, k):
    while True:
        A = []
        for i in range(M):
            if i < k:
                A.append(rng.randrange(1, 1 << k))
            else:
                A.append(rng.randrange(1, 1 << M))
        if rank(A) == M:
            break
    S = [0] * M
    for i in range(M):
        for j in range(i, M):
            if i < k and j < k:
                continue
            if rng.random() < 0.5:
                S[i] |= 1 << j; S[j] |= 1 << i
    return mul(unip(S), levi(A))


_cache = {}


def fast_find_P(rng, X, Y, tries=None):
    sx = _cache.get(X)
    if sx is None:
        sx = _cache[X] = s_of(X)
    sy = _cache.get(Y)
    if sy is None:
        sy = _cache[Y] = s_of(Y)
    if sx[1] != sy[1]:
        return None
    h = rand_stab(rng, sx[1])
    g = mul(inv(sy[0]), mul(h, sx[0]))
    return g


if __name__ == "__main__":
    from sp2 import rand_P, D, neighbors
    rng = random.Random(3)
    import time
    t = time.time()
    for _ in range(300):
        g1 = rand_P(rng); g2 = rand_P(rng)
        X = image(g1, image(D, image(rand_P(rng), image(D, image(rand_P(rng), image(D, L0))))))
        Y = image(g2, X)
        g = fast_find_P(rng, X, Y)
        assert image(g, X) == Y and image(g, L0) == L0
    print("fast_find_P ok, %.2f ms each" % ((time.time() - t) / 300 * 1000))
