"""v-xw3-extract (skeptic, 14 Oct 2026): my own implementation of the XW block-3 letters and paths, written from the
definitions in xw3lib.py's docstring and path_ext.py's docstring, without importing either.

State: complex vector over addresses x = sum_j x_j 8^j (column j = 3 bits, bit r = row r), f columns, S inputs
(shape (8^f, S)).  A letter is a pair (perm on addresses, Z4 phase as a function of the INPUT address):
(g v)[perm(x)] = i^phase(x) v[x].  Words '.'-joined act LEFT TO RIGHT (first letter first).  Applied directly to
vectors letter by letter (no monomial composition), inverses letter by letter in reverse.
D_z: on every column, v <- a v + b v[x ^ (z << 3j)], a = (1+i)/2, b = (1-i)/2.  Move (z, C) = C D_z C^-1.
"""
import numpy as np

A_ = (1 + 1j) / 2
B_ = (1 - 1j) / 2
IP = np.array([1, 1j, -1, -1j])


def bits(c, r):
    return (c >> r) & 1


def par(a):
    a = np.asarray(a); p = np.zeros_like(a)
    while np.any(a):
        p ^= a & 1; a = a >> 1
    return p


class F:
    """addresses at one f"""

    def __init__(self, f):
        self.f = f; self.n = 8 ** f
        x = np.arange(self.n, dtype=np.int64); self.x = x
        self.c = [(x >> (3 * j)) & 7 for j in range(f)]      # column values

    def pack(self, cols):
        out = np.zeros(self.n, np.int64)
        for j, cj in enumerate(cols): out |= (cj.astype(np.int64) & 7) << (3 * j)
        return out

    # ---- letters: return (perm, phase) with perm, phase arrays over input address x
    def col_map(self, fun_perm=None, fun_ph=None, js=None):
        js = range(self.f) if js is None else js
        cols = [cj.copy() for cj in self.c]; ph = np.zeros(self.n, np.int64)
        for j in js:
            a = self.c[j]
            if fun_ph is not None: ph = ph + fun_ph(a)
            if fun_perm is not None: cols[j] = fun_perm(a)
        return self.pack(cols), ph % 4

    def letter(self, nm):
        f = self.f
        if nm == "1": return self.x.copy(), np.zeros(self.n, np.int64)
        base = nm; js = None
        if nm[0] in "lf" and len(nm) > 1 and nm[1] in "XSZCPT":
            js = [f - 1] if nm[0] == "l" else [0]; base = nm[1:]
        b0 = lambda a: bits(a, 0); b1 = lambda a: bits(a, 1); b2 = lambda a: bits(a, 2)
        bb = [b0, b1, b2]
        if base in ("P01", "P12", "P02"):
            r, s = int(base[1]), int(base[2])
            def pm(a, r=r, s=s):
                ar, as_ = bits(a, r), bits(a, s)
                return (a & ~((1 << r) | (1 << s))) | (ar << s) | (as_ << r)
            return self.col_map(pm, None, js)
        if base[0] == "C" and len(base) == 3 and base[1:].isdigit():
            r, s = int(base[1]), int(base[2])
            return self.col_map(lambda a, r=r, s=s: a ^ (bits(a, r) << s), None, js)
        if base[0] == "X" and base[1:].isdigit():
            z = int(base[1:]); return self.col_map(lambda a, z=z: a ^ z, None, js)
        if base[0] == "Z" and base[1:].isdigit():
            z = int(base[1:]); return self.col_map(None, lambda a, z=z: 2 * par(a & z), js)
        if base[0] == "S" and base[1:].isdigit():
            z = int(base[1:]); return self.col_map(None, lambda a, z=z: par(a & z), js)
        if base.startswith("CZ"):
            r, s = int(base[2]), int(base[3]); return self.col_map(None, lambda a, r=r, s=s: 2 * (bits(a, r) & bits(a, s)), js)
        if base.startswith("CS"):
            r, s = int(base[2]), int(base[3]); k = 3 if base.endswith("i") else 1
            return self.col_map(None, lambda a, r=r, s=s, k=k: k * (bits(a, r) & bits(a, s)), js)
        if base in ("T0", "T1", "T2"):
            t = int(base[1]); o = [u for u in range(3) if u != t]
            return self.col_map(lambda a, t=t, o=o: a ^ ((bits(a, o[0]) & bits(a, o[1])) << t), None, js)
        if nm in ("A", "Ai", "B", "Bi") or (nm[0] == "a" and len(nm) in (2, 3)) or (nm[0] == "E" and len(nm) == 3):
            return self.chain_letter(nm)
        if nm.startswith("Pr"):
            r = int(nm[2]); p = np.zeros(self.n, np.int64)
            for j in range(f): p ^= bits(self.c[j], r)
            return self.x.copy(), p % 4
        raise KeyError(nm)

    def chain_letter(self, nm):
        """GF(2)-linear row maps along the columns.  A: x_{r,j} ^= x_{r,j-1} (old values), Ai: prefix XOR along j;
        B / Bi the mirror (x_{r,j} ^= x_{r,j+1}, suffix XOR); a<r>, a<r>i: row r only; E<r><s>: x_{r,j} ^= x_{s,j-1}."""
        f = self.f
        X = [[bits(self.c[j], r) for j in range(f)] for r in range(3)]   # X[r][j] old bits
        Y = [[X[r][j].copy() for j in range(f)] for r in range(3)]
        if nm in ("A", "B") or (nm[0] == "a" and len(nm) == 2):
            rows = range(3) if nm in ("A", "B") else [int(nm[1])]
            for r in rows:
                for j in range(f):
                    src = j - 1 if nm != "B" else j + 1
                    if 0 <= src < f: Y[r][j] = X[r][j] ^ X[r][src]
        elif nm in ("Ai", "Bi") or (nm[0] == "a" and len(nm) == 3):
            rows = range(3) if nm in ("Ai", "Bi") else [int(nm[1])]
            order = list(range(f)) if nm != "Bi" else list(range(f - 1, -1, -1))
            for r in rows:
                acc = np.zeros(self.n, np.int64)
                for j in order:
                    acc = acc ^ X[r][j]; Y[r][j] = acc.copy()
        elif nm[0] == "E":
            r, s = int(nm[1]), int(nm[2])
            for j in range(1, f): Y[r][j] = X[r][j] ^ X[s][j - 1]
        cols = [Y[0][j] | (Y[1][j] << 1) | (Y[2][j] << 2) for j in range(f)]
        return self.pack(cols), np.zeros(self.n, np.int64)


class P:
    """path machinery at one f with S inputs"""

    def __init__(self, f, S, seed):
        self.F = F(f); self.f = f; self.n = self.F.n
        rng = np.random.default_rng(seed)
        self.v = (rng.standard_normal((self.n, S)) + 1j * rng.standard_normal((self.n, S))) / np.sqrt(2)
        self.lc = {}

    def L(self, nm):
        if nm not in self.lc:
            pm, ph = self.F.letter(nm)
            assert np.array_equal(np.sort(pm), self.F.x), nm
            self.lc[nm] = (pm, ph)
        return self.lc[nm]

    def let(self, nm, V, inv=False):
        pm, ph = self.L(nm)
        if not inv:
            Y = np.empty_like(V); Y[pm] = IP[ph][:, None] * V; return Y
        return IP[(-ph) % 4][:, None] * V[pm]

    def word(self, w, V, inv=False):
        ls = [x for x in w.split(".") if x and x != "1"]
        if inv: ls = ls[::-1]
        for nm in ls: V = self.let(nm, V, inv)
        return V

    def D(self, z, V):
        for j in range(self.f):
            V = A_ * V + B_ * V[self.F.x ^ (z << (3 * j))]
        return V

    def move(self, z, C, V):
        if C in ("", "1"): return self.D(z, V)
        return self.word(C, self.D(z, self.word(C, V, inv=True)))


def parse_moves(s):
    out = []
    for m in s.split(";"):
        p = m.split(","); out.append((int(p[0]), p[1] if len(p) > 1 else "1"))
    return out


def top_labels(n, R, W):
    """same order as a nested loop over mix before move 1, ..., mix before move n-1, readout (innermost)"""
    import itertools
    return [("1",) + t for t in itertools.product(*W, R)]


def path_vec(Pp, mv, sub, mixes, front):
    """front . M_{sub[-1]} . mix_{last} ... M_{sub[0]} . mix_{first} . v"""
    V = Pp.v
    for i, k in zip(sub, mixes):
        V = Pp.move(mv[i][0], mv[i][1], Pp.word(k, V))
    return Pp.word(front, V)


def lower_list(n, R, W, B):
    """every sub-path with at least one move dropped: (sub, mixes, front)"""
    import itertools
    out = []
    for k in range(n):
        for sub in itertools.combinations(range(n), k):
            front = R if (k > 0 and sub[-1] == n - 1) else B
            mixl = [["1"] if i == 0 else W[i - 1] for i in sub]
            for mix in itertools.product(*mixl):
                for g in front: out.append((sub, mix, g))
    return out
