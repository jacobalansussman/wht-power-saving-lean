"""xw2lib.py -- mw2-xwlayers (14 Oct 2026).  The everything-free game XW at block m (default 2), self-contained,
modelled on ../ms-xw3/xw3lib.py (block 3) with the block size as a parameter.

Address: column j = m bits (row r = bit r of the column); x = sum_j x_j 2^(m j); f digit columns, N = 2^(m f).
Monomial letter g: (g v)[sig[x]] = i^ph[x] v[x].  Words are '.'-joined and act LEFT TO RIGHT (first letter first).
Paid letters inside words: 'D1', 'D2', 'D3' (= dir(z) on every column, z = 1, 2, 3 at block 2), 'Di1'.. inverses.
Free letters (block 2):
  column-uniform: 1, P (row swap), C01 (x1 ^= x0), C10, X1 X2 X3, Z1 Z2 Z3, S1 S2 S3 (S_z = i^(z.x)), Om (= i per
     column), CZ, CS (= i^(x0 x1) per column), CSi, Neg (-1)
  chain maps: A / Ai (both rows ^= row << 1 / prefix XOR), a0 a1 (one row) and a0i a1i, E01 (row0 ^= row1 << 1), E10,
     B / Bi (mirror: row ^= row >> 1 / suffix XOR)
  parity phases: Pr0, Pr1 = i^(parity of row r)
  boundary (last column 'l', first column 'f'): lX1 lX2 lZ1 lZ2 lS1 lS2 lCZ lCS lC01, and f... alike
"""
import numpy as np

A_C, B_C = (1 + 1j) / 2, (1 - 1j) / 2
PH = np.array([1, 1j, -1, -1j])


class Mono:
    __slots__ = ("sig", "ph")

    def __init__(self, sig, ph):
        self.sig = sig; self.ph = ph

    def then(self, other):                       # other o self (self acts first)
        return Mono(other.sig[self.sig], (self.ph + other.ph[self.sig]) % 4)

    def inv(self):
        sig = np.empty_like(self.sig); sig[self.sig] = np.arange(self.sig.size)
        ph = np.empty_like(self.ph); ph[self.sig] = (-self.ph) % 4
        return Mono(sig, ph)


def col_tables(m):
    """per-column tables name -> (perm or None, Z4 phase or None) on F_2^m"""
    T = {}
    q = 1 << m
    a = np.arange(q)
    b = [(a >> r) & 1 for r in range(m)]
    par = lambda v: np.array([bin(int(t)).count("1") & 1 for t in v])
    for z in range(1, q):
        T["X%d" % z] = (a ^ z, None)
        T["Z%d" % z] = (None, 2 * par(a & z))
        T["S%d" % z] = (None, par(a & z))
    T["Om"] = (None, np.ones(q, np.int64))
    T["Neg"] = (None, 2 * np.ones(q, np.int64))
    if m == 2:
        T["P"] = (((a & 1) << 1) | ((a >> 1) & 1), None)
        T["C01"] = (a ^ (b[0] << 1), None)
        T["C10"] = (a ^ (b[1] << 0), None)
        T["CZ"] = (None, 2 * (b[0] & b[1]))
        T["CS"] = (None, (b[0] & b[1]))
        T["CSi"] = (None, 3 * (b[0] & b[1]))
    return T


class Lett:
    def __init__(self, f, m=2):
        self.f = f; self.m = m; self.q = 1 << m; self.n = self.q ** f
        self.x = np.arange(self.n, dtype=np.int64)
        self.cols = np.stack([(self.x >> (m * j)) & (self.q - 1) for j in range(f)], axis=1)

    def ident(self):
        return Mono(self.x.copy(), np.zeros(self.n, np.int64))

    def fromcols(self, c):
        out = np.zeros(self.n, np.int64)
        for j in range(self.f):
            out |= c[:, j].astype(np.int64) << (self.m * j)
        return out

    def cu(self, perm=None, phase=None, cols_sel=None):
        js = range(self.f) if cols_sel is None else cols_sel
        c = self.cols.copy(); ph = np.zeros(self.n, np.int64)
        for j in js:
            if phase is not None: ph += np.asarray(phase)[self.cols[:, j]]
            if perm is not None: c[:, j] = np.asarray(perm)[self.cols[:, j]]
        return Mono(self.fromcols(c), ph % 4)

    def chain(self, r, s, inverse=False, mirror=False):
        f = self.f; c = self.cols.copy()
        order = list(range(f)) if not mirror else list(range(f - 1, -1, -1))
        if not inverse:
            src = (self.cols >> s) & 1
            for t in range(1, f):
                j, jp = order[t], order[t - 1]
                c[:, j] ^= src[:, jp] << r
        else:
            assert r == s
            acc = np.zeros(self.n, np.int64)
            for t in range(f):
                j = order[t]
                acc ^= (self.cols[:, j] >> r) & 1
                c[:, j] = (self.cols[:, j] & ~(1 << r)) | (acc << r)
        return Mono(self.fromcols(c), np.zeros(self.n, np.int64))

    def parity_phase(self, r):
        return Mono(self.x.copy(), (((self.cols >> r) & 1).sum(axis=1) % 2).astype(np.int64))


class Game:
    """letters at several f; words with paid letters applied to vectors (never dense matrices)."""

    def __init__(self, FS, m=2):
        self.FS = list(FS); self.m = m
        self.L = {f: Lett(f, m) for f in self.FS}
        self.TAB = col_tables(m)
        self.cache = {}

    def letter(self, nm, f):
        T = self.L[f]; TAB = self.TAB
        if nm == "1": return T.ident()
        if nm in TAB:
            pm, ph = TAB[nm]; return T.cu(pm, ph)
        if nm[0] in "lf" and nm[1:] in TAB:
            j = T.f - 1 if nm[0] == "l" else 0
            pm, ph = TAB[nm[1:]]; return T.cu(pm, ph, cols_sel=[j])
        if nm in ("A", "Ai", "B", "Bi"):
            g = T.ident(); inv = nm.endswith("i"); mir = nm[0] == "B"
            for r in range(self.m):
                g = g.then(T.chain(r, r, inverse=inv, mirror=mir))
            return g
        if nm[0] in "ab" and len(nm) in (2, 3) and nm[1].isdigit():
            return T.chain(int(nm[1]), int(nm[1]), inverse=nm.endswith("i"), mirror=nm[0] == "b")
        if nm[0] == "E" and len(nm) == 3:
            return T.chain(int(nm[1]), int(nm[2]))
        if nm.startswith("Pr"):
            return T.parity_phase(int(nm[2]))
        raise KeyError(nm)

    def mono(self, nm):
        """free word (no paid letters) -> {f: Mono}"""
        if nm not in self.cache:
            out = {}
            for f in self.FS:
                g = self.L[f].ident()
                for t in [x for x in nm.split(".") if x and x != "1"]:
                    if t.endswith("^-1"):
                        g = g.then(self.letter(t[:-3], f).inv())
                    else:
                        g = g.then(self.letter(t, f))
                out[f] = g
            self.cache[nm] = out
        return self.cache[nm]

    @staticmethod
    def apply_mono(g, V):
        Y = np.empty_like(V)
        Y[g.sig] = PH[g.ph].reshape((-1,) + (1,) * (V.ndim - 1)) * V
        return Y

    def apply_D(self, z, V, f, inverse=False):
        a, b = (B_C, A_C) if inverse else (A_C, B_C)
        idx = self.L[f].x
        for j in range(f):
            V = a * V + b * V[idx ^ (z << (self.m * j))]
        return V

    def apply_word(self, w, V, f):
        """apply word w (letters left to right, paid letters 'D<z>' / 'Di<z>') to V (n_f x S)"""
        toks = [t for t in w.split(".") if t and t != "1"]
        buf = []
        def flush(V):
            if buf:
                V = self.apply_mono(self.mono(".".join(buf))[f], V); buf.clear()
            return V
        for t in toks:
            if t[0] == "D" and t[1:].isdigit():
                V = flush(V); V = self.apply_D(int(t[1:]), V, f)
            elif t[:2] == "Di" and t[2:].isdigit():
                V = flush(V); V = self.apply_D(int(t[2:]), V, f, inverse=True)
            else:
                buf.append(t)
        return flush(V)

    def K(self, V, f):
        for z in [1 << r for r in range(self.m)]:
            V = self.apply_D(z, V, f)
        return V
