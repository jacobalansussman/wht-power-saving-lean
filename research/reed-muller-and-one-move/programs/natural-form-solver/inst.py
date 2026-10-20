"""inst.py (rm2-solve): instances of the multi-target commutative slice (Theorem B form of the no-phase Singer game).

A program at block m: W arrays, array l owed [T_l] (T_l symmetric over F_2), moves [z_i z_i^T] in time order
i = 0..M-1, scalar mixing free. Exact order-2 system (ex2-families Lemma 1 with a target per array, rm-theory 8.5):
for every binary quadratic form Q (Q(0) = 0), with t_i = Q(z_i) and p_l = <Q, T_l>,
    F(Q) := A T (I + H T)^-1 B^T = diag(p_1..p_W),    T = diag(t), H strictly upper triangular.
Unknowns A, B (W x M) and H (M(M-1)/2): 2WM + M(M-1)/2 complex numbers.

Coordinates: a self-dual basis f_1..f_m of GF(2^m) (rm-theory gf.py), so the dot product is the trace form and the
canonical Singer targets are T_k = sigma^k sigma^kT (frame alpha^k f_j), k in Z_n.
"""
import itertools, random
import numpy as np
from gf import GF


def mons(m):
    return [(i, i) for i in range(m)] + [(i, j) for i in range(m) for j in range(i + 1, m)]


class Field:
    def __init__(self, m, seed=1):
        self.m = m
        F = self.F = GF(m)
        self.fb = F.self_dual_basis(seed)
        self.co = {x: F.coords(x, self.fb) for x in range(F.q)}
        self.toint = {x: sum(c << j for j, c in enumerate(self.co[x])) for x in range(F.q)}

    def frame(self, k):
        """frame alpha^k f_j as coordinate ints"""
        return [self.toint[self.F.mul(self.F.alpha(k), f)] for f in self.fb]

    def sym_of_frame(self, frame):
        m = self.m
        S = [[0] * m for _ in range(m)]
        for v in frame:
            for i in range(m):
                for j in range(m):
                    S[i][j] ^= ((v >> i) & 1) & ((v >> j) & 1)
        return S

    def target(self, k):
        return self.sym_of_frame(self.frame(k))

    def log_of_int(self, z):
        inv = {v: x for x, v in self.toint.items()}
        return self.F.log[inv[z]]


def sym_identity(m):
    return [[int(i == j) for j in range(m)] for i in range(m)]


class Instance:
    """targets: list of m x m symmetric 0/1 matrices; dirs: list of coordinate ints (time order)."""
    def __init__(self, m, targets, dirs, label=""):
        self.m, self.T, self.dirs, self.label = m, targets, list(dirs), label
        self.W, self.M = len(targets), len(dirs)
        ML = mons(m); P = len(ML)
        coeffs = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)   # all 2^P forms
        zb = np.array([[(z >> j) & 1 for j in range(m)] for z in self.dirs], dtype=np.int64)
        mon = np.stack([zb[:, i] * zb[:, j] for (i, j) in ML], axis=1)                 # M x P
        self.t_all = (coeffs @ mon.T) % 2                                                # Nf x M
        tv = np.array([[S[i][j] for (i, j) in ML] for S in targets], dtype=np.int64)    # W x P
        self.p_all = (coeffs @ tv.T) % 2                                                 # Nf x W
        D = self.t_all.sum(1)
        bad = (D == 0) & (self.p_all.sum(1) > 0)
        self.span_ok = not bad.any()
        keep = D > 0
        self.t = self.t_all[keep].astype(float); self.p = self.p_all[keep].astype(float)
        self.coeffs = coeffs[keep]
        self.Nf = self.t.shape[0]

    def exact_tests(self):
        D = self.t.sum(1); nodd = self.p.sum(1)
        return dict(span_ok=self.span_ok, violations=int((D < nodd).sum()),
                    tight=int(((D == nodd) & (D > 0)).sum()))


def canonical(m, ks, dirs="all", order="lex", seed=0, fseed=1):
    """canonical Singer targets T_k for k in ks; dirs 'all' (one move per direction) or 'frames' (plain schedule)."""
    fd = Field(m, fseed)
    T = [fd.target(k) for k in ks]
    if dirs == "all":
        d = list(range(1, 1 << m))
    elif dirs == "frames":
        d = [z for k in ks for z in fd.frame(k)]
    elif dirs == "Iframes":
        d = [z for _ in ks for z in fd.frame(0)]
    else:
        d = list(dirs)
    d = apply_order(fd, d, order, seed)
    return Instance(m, T, d, label=f"m{m} ks={ks} dirs={dirs} order={order}"), fd


def slice_I(m, W, dirs="all", order="lex", seed=0, fseed=1):
    fd = Field(m, fseed)
    T = [sym_identity(m) for _ in range(W)]
    if dirs == "all":
        d = list(range(1, 1 << m))
    elif dirs == "frames":
        d = [1 << j for _ in range(W) for j in range(m)]
    else:
        d = list(dirs)
    d = apply_order(fd, d, order, seed)
    return Instance(m, T, d, label=f"m{m} slice I x{W} dirs={dirs} order={order}"), fd


def apply_order(fd, d, order, seed):
    d = list(d)
    if order == "lex":
        d = sorted(d)
    elif order == "keep":
        pass
    elif order == "singer":
        d = sorted(d, key=fd.log_of_int)
    elif order.startswith("rand"):
        s = int(order[4:]) if order[4:] else seed
        random.Random(s).shuffle(d)
    else:
        raise ValueError(order)
    return d
