"""keptlib.py -- xw-border (14 Oct 2026).  Block 2, XW, with the paid letters replaced by FIBREWISE SUB-CALLS.

A paid move D1 = dir(e1)^(x)f acts on the f row-0 bits (bit 0 of every column) and is the identity on the row-1 bits.
The recursion computes it by a sub-call on those f bits, applied to every row-1 fibre: its value is Y (x) I for some
operator Y on C^(2^f) (the exact Q = sqrt(X)^(x)f, or any approximation Q + eps E1 + ...).  Likewise D2 -> I (x) Y',
and D3 = C01 D1 C01 (C01: x1 ^= x0 in every column) -> C01 (Y'' (x) I) C01.
Vectors: length 4^f, index x = sum_j (x0_j + 2 x1_j) 4^j (as in xw2lib).  Reshaped to M[i0, i1] with
i0 = sum_j x0_j 2^j, i1 = sum_j x1_j 2^j, so (Y (x) I) acts as M -> Y M and (I (x) Y') as M -> M Y'^T.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import xw2lib

A_C, B_C = (1 + 1j) / 2, (1 - 1j) / 2
SQX = np.array([[A_C, B_C], [B_C, A_C]])


def kron_pow(g, f):
    out = np.ones((1, 1), complex)
    for _ in range(f):
        out = np.kron(g, out)          # bit j of the index = j-th factor from the RIGHT (little endian): kron(g, out)
    return out


class Fib:
    def __init__(self, f):
        self.f = f; self.n = 4 ** f
        x = np.arange(self.n)
        i0 = np.zeros(self.n, np.int64); i1 = np.zeros(self.n, np.int64)
        for j in range(f):
            i0 |= ((x >> (2 * j)) & 1) << j
            i1 |= ((x >> (2 * j + 1)) & 1) << j
        self.i0, self.i1 = i0, i1
        self.G = xw2lib.Game([f])
        self.Q = kron_pow(SQX, f)

    def to_mat(self, v):                    # v: (n, S) -> (2^f, 2^f, S)
        M = np.zeros((2 ** self.f, 2 ** self.f) + v.shape[1:], complex)
        M[self.i0, self.i1] = v
        return M

    def to_vec(self, M):
        return M[self.i0, self.i1]

    def D(self, z, v, Y=None):
        """paid letter z applied to v with sub-call value Y (None = exact Q)."""
        if Y is None:
            Y = self.Q
        if z == 1:
            return self.to_vec(np.einsum('ab,bc...->ac...', Y, self.to_mat(v)))
        if z == 2:
            return self.to_vec(np.einsum('cb,ab...->ac...', Y, self.to_mat(v)))
        if z == 3:
            c = self.free("C01", v)
            return self.free("C01", self.D(1, c, Y))
        raise ValueError(z)

    def free(self, w, v):
        return self.G.apply_word(w, v, self.f)

    def rand_Y(self, rng, scale=1.0):
        d = 2 ** self.f
        return scale * (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))) / np.sqrt(2 * d)

    def eq_last_exact_Y(self, rng):
        """a sub-call value that (i) commutes with X^(x)f and (ii) is EXACT sqrt(X) on the last bit:
        Y = Y'' (x) sqrt(X) with Y'' commuting with X^(x)(f-1)."""
        f = self.f
        d = 2 ** (f - 1)
        R = (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))) / np.sqrt(2 * d)
        Xf = kron_pow(np.array([[0, 1], [1, 0]], complex), f - 1)
        R = R + Xf @ R @ Xf
        return np.kron(SQX, R)                 # last bit = most significant = leftmost kron factor
