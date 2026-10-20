"""xw3lib.py -- ms-xw3 (moonshot round, 14 Oct 2026).  The everything-free game XW at BLOCK 3, self-contained.

Address of column j = 3 bits (row r = bit r); array index x = sum_j x_j 8^j; f digit columns.
Monomial letter g: (g v)[sig[x]] = i^ph[x] v[x]  (same convention as wave17 xclib.Mono).
Paid move D_z = dir(z) on every column, dir(z) = a I + b X^z, a = (1+i)/2, b = (1-i)/2.  K = D_1 D_2 D_4.
Spread move: C D_z C^-1 for a free word C.

Letters (names; words are '.'-joined, applied LEFT TO RIGHT, first letter acts first):
 column-uniform (every column):
   1, P01 P12 P02 (row swaps), Crs (x_s ^= x_r, r != s), X<z> Z<z> S<z> (z = 1..7; S<z> = i^[z.x]), Om (= i per column),
   CZrs = (-1)^(x_r x_s), CSrs = i^(x_r x_s), CSrsi (inverse), CCZ = (-1)^(x0 x1 x2),
   Tt  = Toffoli with target t (x_t ^= product of the other two bits)                     [non-Clifford]
 chain maps: A / Ai (every row ^= row << 1 / prefix XOR), a<r> / a<r>i (one row), E<r><s> (row r ^= row s << 1),
             B / Bi (mirror: every row x_j ^= x_(j+1) / suffix XOR)
 parity phases: Pr<r> = i^(parity of row r) = a + b Zstr_r
 boundary (last column 'l', first column 'f'): lX<z>, lS<z>, lZ<z>, lCS01, lCS02, lCS12, lT2, lT0, lT1, and f... alike
"""
import numpy as np

M3 = 3
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


def _par(a):
    a = np.asarray(a); c = np.zeros_like(a)
    for k in range(8):
        c ^= (a >> k) & 1
    return c


class Lett:
    """letters at one f"""

    def __init__(self, f):
        self.f = f; self.n = 8 ** f
        self.x = np.arange(self.n, dtype=np.int64)
        self.cols = np.stack([(self.x >> (3 * j)) & 7 for j in range(f)], axis=1)

    def ident(self):
        return Mono(self.x.copy(), np.zeros(self.n, np.int64))

    def fromcols(self, c):
        out = np.zeros(self.n, np.int64)
        for j in range(self.f):
            out |= c[:, j].astype(np.int64) << (3 * j)
        return out

    def cu(self, perm=None, phase=None, cols_sel=None):
        """per-column (all columns, or only cols_sel) permutation table perm[a] and Z4 phase table phase[a]
        (phase applied as a function of the INPUT address)."""
        js = range(self.f) if cols_sel is None else cols_sel
        c = self.cols.copy(); ph = np.zeros(self.n, np.int64)
        for j in js:
            if phase is not None: ph += np.asarray(phase)[self.cols[:, j]]
            if perm is not None: c[:, j] = np.asarray(perm)[self.cols[:, j]]
        return Mono(self.fromcols(c), ph % 4)

    def chain(self, r, s, inverse=False, mirror=False):
        """row r ^= row s << 1 (x_{r,j} += x_{s,j-1}); r == s & inverse: prefix XOR. mirror: look the other way."""
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


def tables():
    """column tables: name -> (perm or None, phase or None)"""
    T = {}
    a = np.arange(8)
    b = [(a >> r) & 1 for r in range(3)]
    def swap(r, s):
        return (a & ~((1 << r) | (1 << s))) | (b[r] << s) | (b[s] << r)
    T["P01"] = (swap(0, 1), None); T["P12"] = (swap(1, 2), None); T["P02"] = (swap(0, 2), None)
    for r in range(3):
        for s in range(3):
            if r != s:
                T["C%d%d" % (r, s)] = (a ^ (b[r] << s), None)                 # x_s ^= x_r
    for z in range(1, 8):
        T["X%d" % z] = (a ^ z, None)
        T["Z%d" % z] = (None, 2 * _par(a & z))
        T["S%d" % z] = (None, _par(a & z))
    T["Om"] = (None, np.ones(8, np.int64))
    for (r, s) in ((0, 1), (0, 2), (1, 2)):
        T["CZ%d%d" % (r, s)] = (None, 2 * (b[r] & b[s]))
        T["CS%d%d" % (r, s)] = (None, (b[r] & b[s]))
        T["CS%d%di" % (r, s)] = (None, 3 * (b[r] & b[s]))
    T["CCZ"] = (None, 2 * (b[0] & b[1] & b[2]))
    T["T2"] = (a ^ ((b[0] & b[1]) << 2), None)
    T["T0"] = (a ^ ((b[1] & b[2]) << 0), None)
    T["T1"] = (a ^ ((b[0] & b[2]) << 1), None)
    return T


TAB = tables()
INV_NAME = {}
for nm in TAB:
    if nm.startswith("S") and nm[1:].isdigit(): INV_NAME[nm] = None   # handled numerically
NONCLIFF = ("CS", "CCZ", "T0", "T1", "T2")


class Game:
    """letters for several f; words cached; complex screening arithmetic."""

    def __init__(self, FS, W=2, S=1, seed=1):
        self.FS = list(FS); self.W = W; self.S = S
        self.L = {f: Lett(f) for f in self.FS}
        self.n = {f: self.L[f].n for f in self.FS}
        rng = np.random.default_rng(seed)
        self.R = [{f: (rng.standard_normal((self.n[f], S)) + 1j * rng.standard_normal((self.n[f], S))) / np.sqrt(2)
                   for f in self.FS} for _ in range(W)]
        self.cache = {}
        self.dcache = {}

    # ---------------------------------------------------------------- letters
    def letter(self, nm, f):
        T = self.L[f]
        if nm == "1": return T.ident()
        if nm in TAB:
            pm, ph = TAB[nm]; return T.cu(pm, ph)
        if nm[0] in "lf" and nm[1:] in TAB:
            j = T.f - 1 if nm[0] == "l" else 0
            pm, ph = TAB[nm[1:]]; return T.cu(pm, ph, cols_sel=[j])
        if nm == "A" or nm == "Ai" or nm == "B" or nm == "Bi":
            g = T.ident(); inv = nm.endswith("i"); mir = nm[0] == "B"
            for r in range(3):
                g = g.then(T.chain(r, r, inverse=inv, mirror=mir))
            return g
        if nm[0] == "a" and len(nm) in (2, 3):
            return T.chain(int(nm[1]), int(nm[1]), inverse=nm.endswith("i"))
        if nm[0] == "b" and len(nm) in (2, 3):
            return T.chain(int(nm[1]), int(nm[1]), inverse=nm.endswith("i"), mirror=True)
        if nm[0] == "E" and len(nm) == 3:
            return T.chain(int(nm[1]), int(nm[2]))
        if nm.startswith("Pr"):
            return T.parity_phase(int(nm[2]))
        if nm == "Neg":
            return Mono(T.x.copy(), np.full(T.n, 2, np.int64))
        raise KeyError(nm)

    def word(self, w):
        if w not in self.cache:
            out = {}
            for f in self.FS:
                g = self.L[f].ident()
                for nm in [x for x in w.split(".") if x]:
                    if nm.endswith("^-1"):
                        g = g.then(self.word(nm[:-3])[f].inv() if False else self._inv_letter(nm[:-3], f))
                    else:
                        g = g.then(self.letter(nm, f))
                out[f] = g
            self.cache[w] = out
        return self.cache[w]

    def _inv_letter(self, nm, f):
        return self.letter(nm, f).inv()

    def winv(self, w):
        key = "INV(" + w + ")"
        if key not in self.cache:
            g = self.word(w)
            self.cache[key] = {f: g[f].inv() for f in self.FS}
        return self.cache[key]

    # ---------------------------------------------------------------- actions on contents dict f -> (n_f, ...)
    def mono(self, g, X):
        out = {}
        for f in self.FS:
            mo = g[f]; V = X[f]; Y = np.empty_like(V)
            ph = PH[mo.ph]
            Y[mo.sig] = ph.reshape((-1,) + (1,) * (V.ndim - 1)) * V
            out[f] = Y
        return out

    def W_(self, w, X):
        if w in ("", "1"): return X
        return self.mono(self.word(w), X)

    def D(self, z, X, inverse=False):
        a, b = (B_C, A_C) if inverse else (A_C, B_C)
        out = {}
        for f in self.FS:
            V = X[f]; idx = self.L[f].x
            for j in range(f):
                V = a * V + b * V[idx ^ (z << (3 * j))]
            out[f] = V
        return out

    def move(self, z, C, X):
        if C in ("", "1"): return self.D(z, X)
        return self.mono(self.word(C), self.D(z, self.mono(self.winv(C), X)))

    def K(self, X):
        for z in (1, 2, 4):
            X = self.D(z, X)
        return X

    def inp(self, l):
        return {f: self.R[l][f].copy() for f in self.FS}


def add(X, Y, c=1.0):
    return {f: X[f] + c * Y[f] for f in X}


def flat(X):
    return np.concatenate([X[f].reshape(-1) for f in sorted(X)])


# ------------------------------------------------------------------------ menus
ENG = ["1", "P01", "P12", "C01", "C12", "C20", "S1", "S2", "S4", "S3", "X1", "X2", "X4", "Z1", "Z2", "Z4", "Om",
       "CZ01", "CZ12"]
CHAIN = ["A", "Ai", "a0", "a0i", "a1", "a1i", "E01", "E12", "B", "Bi"]
BDRY = ["lS1", "lS2", "lX1", "lX4", "fS1", "fX1", "Pr0", "Pr1", "Pr2"]
NC = ["T2", "T0", "CS01", "CS12", "CS01i", "CCZ", "lT2", "lCS01", "fT2"]
NCCH = ["A.T2.Ai", "Ai.T2.A", "A.CS01.Ai", "Ai.CS01.A", "E01.T2.E01", "a0.T2.a0i", "B.CS01.Bi", "T2.A.T2"]
# matched Clifford stand-ins (same places, Clifford letters): CNOT for Toffoli, CZ for CS
CL = ["C02", "C10", "CZ01.S1", "CZ12.S2", "CZ01.S1.S1.S1", "CZ02", "lC02", "lCZ01", "fC02"]
CLCH = ["A.C02.Ai", "Ai.C02.A", "A.CZ01.Ai", "Ai.CZ01.A", "E01.C02.E01", "a0.C02.a0i", "B.CZ01.Bi", "C02.A.C02"]


LONGORD = ["A.T2.B", "A.B", "A.C01.B", "a0.T2.b0", "A.T0.B", "A.CS01.B", "B.T2.A", "a0.b0", "A.T2.B.T0"]


def menu(kind):
    if kind.endswith("+"):                     # e.g. 'xw+' = xw plus 60 two-word products
        return menu_ext(kind[:-1])
    if kind == "eng": return ENG + ["C02", "C10", "CZ02"]
    if kind == "longord":                       # free words of exponential order (two-sided chains, with/without Toffoli)
        return LONGORD
    if kind == "xwg": return ENG + CHAIN + BDRY + NC + NCCH + LONGORD + [w + ".X1" for w in LONGORD] + [w + ".Z1" for w in LONGORD]
    if kind == "xq": return ENG + CHAIN + BDRY
    if kind == "xw": return ENG + CHAIN + BDRY + NC + NCCH
    if kind == "xqm": return ENG + CHAIN + BDRY + CL + CLCH
    if kind == "wd": return ENG + NC[:6] + ["C02", "C10", "CZ02"]          # column-uniform Toffoli/CS, no chain
    raise ValueError(kind)


def menu_ext(kind, n_extra=60, seed=7):
    """menu + n_extra random two-word products of its words (a larger slice of the free algebra)."""
    import random
    base = menu(kind)
    rng = random.Random(seed)
    ext = set()
    while len(ext) < n_extra:
        a, b = rng.choice(base[1:]), rng.choice(base[1:])
        if a != b: ext.add(a + "." + b)
    return base + sorted(ext)
