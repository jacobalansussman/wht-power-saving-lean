"""o3lib.py -- o3-build (moonshot o3, 14 Oct 2026): exact word engine for block-4 templates with order-3 letters.

Game (register form, as tb4-beam bm.py, WIDENED letters = every monomial matrix with entries in {1,i,-1,-i}; the
order-3 letters used here are address maps, so they are also engine letters):
  generators start as the inputs x_l; a move makes u = D1 b c, c = sum_j coef_j * a_j * G_j (a_j free, coef rational);
  word correct when K x_l lies in the Q[B']-span of the generators (span test), or, for a PRIVATE END, when
  K x_l in span + E_l * span for one paid element E_l (bm.Game.free_last_test).  2-in-7 with the forced private end
  = a 5-move prefix whose two arrays both pass free_last_test (E1 = E2 is the halving end K x_{1,2} = E(c1 +- c2)).
Program record (JSON-able): list of moves {"b": [sigma, q], "cs": [[num, den, [sigma, q] or None, j], ...]}.
Letters as monomials (sigma, q): e_y -> i^q[y] e_sigma[y].
"""
import sys, os, itertools, random, json, time
PYLIB = ''
if os.path.isdir(PYLIB):
    sys.path.insert(0, PYLIB)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from fractions import Fraction
from bm import Game

M_ = 4
N_ = 1 << M_

def mono_tuple(G, M):
    s, q = G.mono_of(M)
    return [int(x) for x in s], [int(x) for x in q]

def perm_bits(perm, m=M_):
    """address map moving bit j to bit perm[j]"""
    return [sum(((x >> j) & 1) << perm[j] for j in range(m)) for x in range(1 << m)]

def dot(a, b):
    return bin(a & b).count('1') & 1

class W4:
    """block-4 helper around bm.Game"""
    def __init__(self, seed=1):
        self.G = G = Game(M_)
        self.rng = random.Random(seed)
        n = N_
        self.ID = G.I
        self.sig = {}       # order-3 address maps: name -> matrix
        for tri in itertools.combinations(range(M_), 3):
            for o in (0, 1):
                a, b, c = tri if o == 0 else (tri[0], tri[2], tri[1])
                perm = list(range(M_)); perm[a] = b; perm[b] = c; perm[c] = a
                self.sig[f"s{a}{b}{c}"] = G.mono(perm_bits(perm), [0] * n)
        self.Z = {v: G.mono(list(range(n)), [2 * dot(v, x) for x in range(n)]) for v in range(n)}
        self.X = {v: G.mono([x ^ v for x in range(n)], [0] * n) for v in range(n)}
        self.D1 = G.D1
        self.Dinv = np.conj(G.D1).T

    # ---------------------------------------------------------------- generators
    def inp(self, l):
        return self.G.inp(l)

    def lm(self, A, gen):
        return self.G.lmul(A, gen)

    def add(self, *pairs):
        return self.G.add(*pairs)

    def move(self, b, content):
        """u = D1 b content"""
        return self.G.lmul(self.G.mul(self.D1, b), content)

    def depth(self, gid):
        return int(round(np.log2(self.G.rowsupp(self.G.MATS[gid]).max())))

    def dist(self, gid):
        return self.G.dist(gid)

    def profile(self, gen):
        """{(l, depth, dist): count}"""
        out = {}
        for (l, g) in gen:
            k = (l, self.depth(g), self.dist(g))
            out[k] = out.get(k, 0) + 1
        return out

    def leads(self, gen, l=None):
        """geodesic terms of max depth (optionally of array l)"""
        ts = [(ll, g) for (ll, g) in gen if (l is None or ll == l) and self.depth(g) + self.dist(g) == M_]
        if not ts:
            return []
        dm = max(self.depth(g) for (_, g) in ts)
        return [(ll, g) for (ll, g) in ts if self.depth(g) == dm]

    def progress(self, gid, want=6):
        """letters b with D1 b g one step nearer K (full halving chain exists)"""
        return self.G.progress_letters(gid, self.rng, want)

    def advances(self, b, gid):
        g = self.G.MATS[gid]
        h = self.G.mul(self.D1, b, g)
        hid = self.G.intern(h)
        return self.depth(hid) == self.depth(gid) + 1 and self.dist(hid) == self.dist(gid) - 1

    def align(self, T, g):
        """b with D1 b g == T exactly (b monomial) or None; T, g matrices"""
        b = self.G.mul(self.Dinv, T, np.conj(g).T)
        return b if self.G.is_monomial(b) else None

    # ---------------------------------------------------------------- end tests
    def end_private(self, gens, cap=4000):
        """forced private end: each array passes free_last_test (K x_l in span + E_l span). returns (ok, [status], sols)"""
        sts = []; sols = []
        for l in range(2):
            st, sol = self.G.free_last_test(gens, l, cap)
            sts.append(st); sols.append(sol)
        ok = all(s in ("yes", "yes0") for s in sts)
        return ok, sts, sols

    def end_span(self, gens, cap=4000):
        sts = []; sols = []
        for l in range(2):
            st, sol = self.G.span_test(gens, l, cap)
            sts.append(st); sols.append(sol)
        return all(s == "yes" for s in sts), sts, sols

# -------------------------------------------------------------------- programs: record, rebuild, verify
def rec_move(G, b, cs):
    """cs: list of (Fraction coef, letter matrix or None, j)"""
    return {"b": mono_tuple(G, b),
            "cs": [[int(Fraction(c).numerator), int(Fraction(c).denominator),
                    (mono_tuple(G, a) if a is not None else None), j] for (c, a, j) in cs]}

def build(w, prog):
    """rebuild generators from a program record (inputs x0, x1 first)"""
    G = w.G
    gens = [G.inp(0), G.inp(1)]
    for mv in prog:
        b = G.mono(*mv["b"])
        c = {}
        for (nu, de, a, j) in mv["cs"]:
            t = gens[j] if a is None else G.lmul(G.mono(*a), gens[j])
            c = G.add((1, c), (Fraction(nu, de), t))
        gens.append(G.lmul(G.mul(w.D1, b), c))
    return gens

def leftover(w, gens, l, sol):
    """formal value of the solution combination minus K x_l: dict (l2, gid) -> coef; must hold only free terms
    (or E-free terms) for a valid solution"""
    G = w.G
    E, terms = (sol if isinstance(sol, tuple) else (None, sol))
    acc = {}
    for ((kind, i, bid), x) in terms:
        pre = G.MATS[bid] if kind == 0 else G.mul(E, G.MATS[bid])
        for (l2, g2), co in gens[i].items():
            t = (l2, G.intern(G.mul(pre, G.MATS[g2])))
            acc[t] = acc.get(t, 0) + x * co
    acc = {t: v for t, v in acc.items() if v != 0}
    acc[(l, G.KID)] = acc.get((l, G.KID), 0) - 1
    return {t: v for t, v in acc.items() if v != 0}

def formal_check(w, gens, sols):
    """the solution combination equals K x_l up to terms g x_l2 with g free (cancelled by free images of the inputs)
    or E^-1 g free (E applied to free images of the inputs: part of the private last move's content)"""
    G = w.G
    for l, sol in enumerate(sols):
        if sol is None:
            return False
        E = sol[0] if isinstance(sol, tuple) else None
        for (l2, g), v in leftover(w, gens, l, sol).items():
            if G.is_free(g):
                continue
            if E is not None and G.is_monomial(G.mul(np.conj(E).T, G.MATS[g])):
                continue
            return False
    return True

def apply_tensor(M, v, f, n=N_):
    """M^(x)f applied to a vector of length n^f"""
    t = v.reshape((n,) * f)
    for ax in range(f):
        t = np.tensordot(M, t, axes=([1], [ax]))
        t = np.moveaxis(t, 0, ax)
    return t.reshape(-1)

def raw_check(w, prog, gens, sols, fs=(1, 2, 3), seed=5):
    """replay the word on random complex arrays at each f; output_l = (solution combination) - (formal free leftover,
    applied to the inputs, i.e. the free images of the inputs that the word adds) must equal K^(x)f x_l.
    Returns {f: max relative error over both arrays}."""
    G = w.G
    res = {}
    rng = np.random.default_rng(seed)
    for f in fs:
        Nf = N_ ** f
        xs = [rng.normal(size=Nf) + 1j * rng.normal(size=Nf) for _ in range(2)]
        vecs = list(xs)
        for mv in prog:
            c = np.zeros(Nf, dtype=complex)
            for (nu, de, a, j) in mv["cs"]:
                v = vecs[j] if a is None else apply_tensor(G.mono(*a), vecs[j], f)
                c = c + (nu / de) * v
            vecs.append(apply_tensor(G.mul(w.D1, G.mono(*mv["b"])), c, f))
        err = 0.0
        for l, sol in enumerate(sols):
            E, terms = (sol if isinstance(sol, tuple) else (None, sol))
            out = np.zeros(Nf, dtype=complex); cE = np.zeros(Nf, dtype=complex)
            for ((kind, i, bid), x) in terms:
                v = apply_tensor(G.MATS[bid], vecs[i], f)
                if kind == 0:
                    out += float(x) * v
                else:
                    cE += float(x) * v
            if E is not None:
                out += apply_tensor(E, cE, f)
            for (l2, g), v in leftover(w, gens, l, sol).items():
                out -= float(v) * apply_tensor(G.MATS[g], xs[l2], f)
            target = apply_tensor(G.K, xs[l], f)
            err = max(err, float(np.abs(out - target).max() / np.abs(target).max()))
        res[f] = err
    return res
