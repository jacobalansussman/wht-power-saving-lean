"""v-w3-blueprint: own minimal model of the block-m game (independent of symsolve/bm).
Elements = 2^m x 2^m complex matrices, exact on a dyadic grid (snapped). D(z) = ((1+i)I + (1-i)X^z)/2, K = prod D(e_j).
Engine letters B: affine address maps y -> A y + c composed with Z4 phases i^(a.y + 2 Q(y)). Widened: any Z4 monomial."""
import numpy as np, itertools, random
SC = 1 << 14
def snap(M):
    R = np.round(M.real * SC) / SC; I = np.round(M.imag * SC) / SC
    assert np.max(np.abs(R - M.real)) < 1e-9 and np.max(np.abs(I - M.imag)) < 1e-9
    return R + 1j * I
class Grp:
    def __init__(s, m):
        s.m = m; s.n = n = 1 << m; s.I = np.eye(n, dtype=complex)
        s.D = {z: snap(((1 + 1j) * s.I + (1 - 1j) * s.X(z)) / 2) for z in range(1, n)}
        K = s.I
        for j in range(m): K = s.mul(K, s.D[1 << j])
        s.K = K; s.Ki = np.conj(K).T
    def X(s, z):
        M = np.zeros((s.n, s.n), dtype=complex)
        for y in range(s.n): M[y ^ z, y] = 1
        return M
    def mul(s, *Ms):
        out = Ms[0]
        for M in Ms[1:]: out = snap(out @ M)
        return out
    def inv(s, M): return np.conj(M).T
    def mono(s, sig, q):
        M = np.zeros((s.n, s.n), dtype=complex)
        for y in range(s.n): M[sig[y], y] = 1j ** (q[y] % 4)
        return M
    def is_mono(s, M):
        nz = np.abs(M) > 1e-9
        if not ((nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()): return False
        return all(min(abs(v - z) for z in (1, 1j, -1, -1j)) < 1e-9 for v in M[nz])
    def is_mono_any(s, M):   # monomial up to a global phase (any unit scalars on a common phase class)
        nz = np.abs(M) > 1e-9
        return (nz.sum(0) == 1).all() and (nz.sum(1) == 1).all()
    def dist_in(s, M):       # input distance of the cell B M
        return int(round(np.log2((np.abs(M) > 1e-9).sum(1).max())))
    def apply_gl(s, A, y):
        v = 0
        for j in range(s.m):
            if (y >> j) & 1: v ^= A[j]
        return v
    def rand_gl(s, rng):
        while True:
            A = [rng.randrange(1, s.n) for _ in range(s.m)]
            span = {0}; ok = True
            for a in A:
                if a in span: ok = False; break
                span |= {x ^ a for x in span}
            if ok: return A
    def phase_q(s, a, C):
        q = []
        for y in range(s.n):
            b = [(y >> j) & 1 for j in range(s.m)]
            v = sum(a[j] * b[j] for j in range(s.m)) + 2 * sum(C[j][k] * b[j] * b[k] for j in range(s.m) for k in range(j + 1, s.m))
            q.append(v % 4)
        return q
    def rand_engine(s, rng):
        A = s.rand_gl(rng); c = rng.randrange(s.n)
        sig = [s.apply_gl(A, y) ^ c for y in range(s.n)]
        a = [rng.randrange(4) for _ in range(s.m)]
        C = [[rng.randrange(2) for _ in range(s.m)] for _ in range(s.m)]
        q = s.phase_q(a, C)
        Ph = s.mono(list(range(s.n)), q); Pm = s.mono(sig, [0] * s.n)
        return s.mul(Ph, Pm) if rng.random() < .5 else s.mul(Pm, Ph)
    def rand_wide(s, rng):
        sig = list(range(s.n)); rng.shuffle(sig)
        return s.mono(sig, [rng.randrange(4) for _ in range(s.n)])
    def addr(s, A):
        return s.mono([s.apply_gl(A, y) for y in range(s.n)], [0] * s.n)
    def hmax_lift(s, A):
        """c = phase * addr(A) with K^-1 c^-1 K monomial (search over Z4 linear + quadratic phases)"""
        Ad = s.addr(A)
        for Cb in itertools.product([0, 1], repeat=s.m * (s.m - 1) // 2):
            C = [[0] * s.m for _ in range(s.m)]; it = iter(Cb)
            for j in range(s.m):
                for k in range(j + 1, s.m): C[j][k] = next(it)
            for a in itertools.product(range(4), repeat=s.m):
                c = s.mul(s.mono(list(range(s.n)), s.phase_q(a, C)), Ad)
                cp = s.mul(s.Ki, s.inv(c), s.K)
                if s.is_mono(cp): return c
        return None
