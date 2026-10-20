"""mini.py -- v-mw2-xwlayers (skeptic), independent block-2 XW mini-game, written from scratch (not using xw2lib/als2).

Column value c = r0 + 2 r1 (two address bits per digit column), flat index i = sum_j c_j 4^j.
Paid: D_z = dir(z)^{(x) f}, dir(z) = a I + b X^z, a = (1+i)/2, b = (1-i)/2.  Target K = (dir(1) dir(2))^{(x) f}.
Free letters (monomial): CS, CSi, C01 (r1 ^= r0), C10, S0, S1 (Z4 phases), X0, Z0, chi0/chi0i/chi1/chi1i (running XOR of
one row across columns and its inverse), P0/P1 (parity phase of a row: 1 if even, i if odd), fCS/lCS (CS on the first /
last column only), A (= chi0 . CS . chi0i, a column-mixing conjugate).
"""
import numpy as np

a = (1 + 1j) / 2
b = (1 - 1j) / 2


class Game:
    def __init__(self, f):
        self.f = f
        self.N = 4 ** f
        i = np.arange(self.N)
        self.c = [(i >> (2 * j)) & 3 for j in range(f)]
        self.r0 = [cj & 1 for cj in self.c]
        self.r1 = [(cj >> 1) & 1 for cj in self.c]
        self._cache = {}

    # --- index helpers
    def idx(self, r0, r1):
        out = np.zeros(self.N, dtype=np.int64)
        for j in range(self.f):
            out |= (r0[j].astype(np.int64) << (2 * j)) | (r1[j].astype(np.int64) << (2 * j + 1))
        return out

    def perm(self, name):
        """returns array p with (L x)[p[i]] = x[i]"""
        r0, r1, f = list(self.r0), list(self.r1), self.f
        if name == "C01":
            r1 = [r1[j] ^ r0[j] for j in range(f)]
        elif name == "C10":
            r0 = [r0[j] ^ r1[j] for j in range(f)]
        elif name == "X0":
            r0 = [r0[j] ^ 1 for j in range(f)]
        elif name in ("chi0", "chi1"):
            row = r0 if name == "chi0" else r1
            acc = np.zeros(self.N, dtype=np.int64); new = []
            for j in range(f):
                acc = acc ^ row[j]; new.append(acc.copy())
            if name == "chi0": r0 = new
            else: r1 = new
        elif name in ("chi0i", "chi1i"):
            row = r0 if name == "chi0i" else r1
            new = [row[0]] + [row[j] ^ row[j - 1] for j in range(1, f)]
            if name == "chi0i": r0 = new
            else: r1 = new
        else:
            raise KeyError(name)
        return self.idx(r0, r1)

    def diag(self, name):
        f = self.f
        if name == "CS":
            e = sum(self.r0[j] * self.r1[j] for j in range(f)); return 1j ** e
        if name == "CSi":
            e = sum(self.r0[j] * self.r1[j] for j in range(f)); return (-1j) ** e
        if name == "fCS":
            return 1j ** (self.r0[0] * self.r1[0])
        if name == "lCS":
            return 1j ** (self.r0[f - 1] * self.r1[f - 1])
        if name == "S0":
            return 1j ** sum(self.r0)
        if name == "S1":
            return 1j ** sum(self.r1)
        if name == "Z0":
            return (-1.0) ** sum(self.r0)
        if name in ("P0", "P1"):
            row = self.r0 if name == "P0" else self.r1
            par = sum(row) & 1
            return np.where(par == 0, 1.0 + 0j, 1j)
        raise KeyError(name)

    def letter(self, name, X):
        """apply one free letter to X of shape (N, S)"""
        if name in ("1", ""):
            return X
        if name == "A":
            return self.letter("chi0", self.letter("CS", self.letter("chi0i", X)))
        key = ("L", name)
        if key not in self._cache:
            try:
                self._cache[key] = ("d", self.diag(name)[:, None])
            except KeyError:
                self._cache[key] = ("p", self.perm(name))
        kind, data = self._cache[key]
        if kind == "d":
            return data * X
        Y = np.empty_like(X); Y[data] = X; return Y

    def word(self, w, X):
        """w = 'L1.L2.L3' means L1 L2 L3 (rightmost applied first)"""
        if w in ("1", ""):
            return X
        for name in reversed(w.split(".")):
            X = self.letter(name, X)
        return X

    def col_op(self, U, X):
        S = X.shape[1]
        T = X.T.reshape((S,) + (4,) * self.f)          # axis k+1 = column f-1-k
        for ax in range(1, self.f + 1):
            T = np.moveaxis(np.tensordot(U, T, axes=([1], [ax])), 0, ax)
        return T.reshape(S, self.N).T

    @staticmethod
    def dir1(z):
        U = a * np.eye(4, dtype=complex)
        for c in range(4):
            U[c ^ z, c] += b
        return U

    def D(self, z, X):
        return self.col_op(self.dir1(z), X)

    def K(self, X):
        return self.col_op(self.dir1(1) @ self.dir1(2), X)
