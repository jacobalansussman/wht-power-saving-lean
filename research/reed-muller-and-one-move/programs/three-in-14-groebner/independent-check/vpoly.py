"""vpoly.py (v-r3-3in14): my own sparse polynomials over Q and my own equation generator.

Model (rm-theory 8.5 / ex2-families Lemma 1, re-derived in NOTE.md section 2): at form Q with move support D
(in time order), the program needs  X_D = A_D (I + G_D)^-1 B_D^T = diag(p(Q)),  where A_D, B_D are W x |D| (column
of move s: a_s, b_s) and G_D is STRICTLY LOWER triangular in time order (G[s][t], s later than t), free.  I build
(I + G)^-1 = sum_k (-G)^k by matrix arithmetic (not by listing chains, as the track did).
A polynomial is a dict: sorted tuple of variable indices (repeats = powers) -> Fraction."""
from fractions import Fraction
import itertools, re


def padd(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + s * c
        if r[m] == 0: del r[m]
    return r


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(sorted(m1 + m2)); r[m] = r.get(m, 0) + c1 * c2
            if r[m] == 0: del r[m]
    return r


def var(i):
    return {(i,): Fraction(1)}


ONE = {(): Fraction(1)}


def matmul(X, Y):
    n, k, m = len(X), len(Y), len(Y[0])
    Z = [[{} for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            acc = {}
            for t in range(k):
                if X[i][t] and Y[t][j]: acc = padd(acc, pmul(X[i][t], Y[t][j]))
            Z[i][j] = acc
    return Z


class VRing:
    def __init__(self, W, U):
        self.W = W; self.U = list(U); self.names = []
        self.a = {}; self.b = {}; self.g = {}
        for l in range(W):
            for i in U: self.a[l, i] = self._new(f"a{l}_{i}")
        for l in range(W):
            for i in U: self.b[l, i] = self._new(f"b{l}_{i}")
        for i in U:
            for j in U:
                if i != j: self.g[i, j] = self._new(f"g{i}_{j}")

    def _new(self, n):
        self.names.append(n); return len(self.names) - 1


def equations(R, order, sups):
    """sups: list of (Dset as list of move indices, parity tuple p_l).  order: list of the moves of U in time."""
    pos = {s: k for k, s in enumerate(order)}; W = R.W; polys = []
    for D, p in sups:
        S = sorted(D, key=lambda s: pos[s]); n = len(S)
        Gm = [[var(R.g[S[i], S[j]]) if i > j else {} for j in range(n)] for i in range(n)]
        negG = [[{m: -c for m, c in e.items()} for e in row] for row in Gm]
        Inv = [[dict(ONE) if i == j else {} for j in range(n)] for i in range(n)]
        Pw = [row[:] for row in Inv]
        for _ in range(n - 1):
            Pw = matmul(Pw, negG)
            Inv = [[padd(Inv[i][j], Pw[i][j]) for j in range(n)] for i in range(n)]
        A = [[var(R.a[l, s]) for s in S] for l in range(W)]
        Bt = [[var(R.b[l, s]) for l in range(W)] for s in S]
        X = matmul(matmul(A, Inv), Bt)
        for l in range(W):
            for l2 in range(W):
                f = dict(X[l][l2])
                if l == l2 and p[l]: f = padd(f, ONE, -1)
                if f: polys.append(f)
    return polys


def to_sing(p):
    out = []
    for m, c in sorted(p.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        mono = "*".join(f"x({i + 1})" for i in m)
        cs = str(abs(c))
        body = (mono if abs(c) == 1 else cs + "*" + mono) if m else cs
        out.append(("-" if c < 0 else "+") + body)
    s = "".join(out)
    return s[1:] if s.startswith("+") else s


def to_msolve(p, ren):
    out = []
    for m, c in sorted(p.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        # powers written as x^k (msolve 0.6.5 rule); here every monomial is multilinear anyway
        from collections import Counter
        cnt = Counter(m)
        mono = "*".join((f"x{ren[i] + 1}" if e == 1 else f"x{ren[i] + 1}^{e}") for i, e in sorted(cnt.items()))
        cs = str(abs(c))
        body = (mono if abs(c) == 1 else cs + "*" + mono) if m else cs
        out.append(("-" if c < 0 else "+") + body)
    s = "".join(out)
    return s[1:] if s.startswith("+") else s


TERM = re.compile(r"([+-]?)([^+-]+)")


def parse_sing(s):
    """parse a Singular polynomial string (short = 0 style: x(1)*x(2)^2, coefficients like 3/4) -> dict."""
    s = s.replace(" ", "").replace("\n", "")
    if s == "0": return {}
    p = {}
    for sign, body in TERM.findall(s):
        c = Fraction(1); vs = []
        for f in body.split("*"):
            mm = re.fullmatch(r"x\((\d+)\)(?:\^(\d+))?", f)
            if mm:
                vs += [int(mm.group(1)) - 1] * int(mm.group(2) or 1)
            else:
                c *= Fraction(f)
        if sign == "-": c = -c
        m = tuple(sorted(vs)); p[m] = p.get(m, 0) + c
        if p[m] == 0: del p[m]
    return p
