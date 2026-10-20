#!/usr/bin/env python3
"""vengine.py (v-fp-capacity): my own exact engine for the engine game at block m, built from scratch
(no code of fp-capacity, sl-combinatorics or v-slack imported).
- group elements = exact 2^m x 2^m matrices (dyadic Gaussian entries; exactness asserted on a 2^-24 grid);
- contents = formal sums {(group element, input): coefficient};
- Sp(2m,2) action gbar by conjugation of Paulis X^a Z^b (label a | b<<m);
- classes: Sp labels gbar*hbar_l (H = Paulis.scalars) or image frames gbar L0 (H = B, raw gauge);
- liveness: fp-capacity Definition 1.1 (class graph, backward pass), implemented term-wise;
- raw arrays: the word run on random complex arrays with the f-fold action g -> g^(x)f.
"""
import numpy as np, itertools, random

class Eng:
    def __init__(s, m):
        s.m = m; s.n = 1 << m; s.mats = {}; s.spc = {}; s.spi = {}
        s.ca = (1 + 1j) / 2; s.cb = (1 - 1j) / 2
        s.I = np.eye(s.n, dtype=complex)
        s.Kmat = s.I
        for j in range(m): s.Kmat = s.dirm(1 << j) @ s.Kmat
        s.IDk = s.key(s.I); s.Kk = s.key(s.Kmat)
    # ---- letters
    def Xm(s, z):
        P = np.zeros((s.n, s.n), dtype=complex)
        for x in range(s.n): P[x ^ z, x] = 1
        return P
    def Zm(s, v): return np.diag([(-1) ** bin(v & x).count("1") for x in range(s.n)]).astype(complex)
    def dirm(s, z): return s.ca * s.I + s.cb * s.Xm(z)
    def addr(s, perm):          # coordinate permutation: bit i -> bit perm[i]
        P = np.zeros((s.n, s.n), dtype=complex)
        for x in range(s.n):
            y = 0
            for i in range(s.m):
                if (x >> i) & 1: y |= 1 << perm[i]
            P[y, x] = 1
        return P
    def phase(s, z, c=0): return np.diag([1j ** ((bin(z & x).count("1") + c) % 2) for x in range(s.n)])
    # ---- group elements
    def key(s, M):
        R = M * (1 << 24); re = np.round(R.real); im = np.round(R.imag)
        assert np.array_equal(re, R.real) and np.array_equal(im, R.imag), "inexact"
        k = (re.astype(np.int64).tobytes(), im.astype(np.int64).tobytes())
        if k not in s.mats: s.mats[k] = M.copy()
        return k
    def decode(s, M):           # M = c X^a Z^b  ->  a | b<<m
        col = M[:, 0]; a = int(np.flatnonzero(np.abs(col) > 1e-12)[0]); c0 = M[a, 0]; b = 0
        for j in range(s.m):
            r = M[a ^ (1 << j), 1 << j] / c0
            assert abs(abs(r) - 1) < 1e-9 and abs(r.imag) < 1e-9
            if r.real < 0: b |= 1 << j
        # verify fully
        P = c0 * s.Xm(a) @ s.Zm(b)
        assert np.allclose(P, M), "not a Pauli"
        return a | (b << s.m)
    def sp_of_mat(s, g):
        gi = g.conj().T
        assert np.allclose(gi @ g, s.I)
        im = []
        for j in range(2 * s.m):
            v = 1 << j; a, b = v & (s.n - 1), v >> s.m
            im.append(s.decode(g @ (s.Xm(a) @ s.Zm(b)) @ gi))
        return tuple(im)
    def sp(s, k):
        if k not in s.spc: s.spc[k] = s.sp_of_mat(s.mats[k])
        return s.spc[k]
    def spinv(s, k):
        if k not in s.spi: s.spi[k] = s.sp_of_mat(s.mats[k].conj().T)
        return s.spi[k]
    @staticmethod
    def apply(S, v):
        out = 0; j = 0
        while v:
            if v & 1: out ^= S[j]
            v >>= 1; j += 1
        return out
    def comp(s, A, B): return tuple(s.apply(A, b) for b in B)    # A o B
    def frame(s, S, base):      # span of S(base vectors)
        F = {0}
        for b in base:
            w = s.apply(S, b); F |= {x ^ w for x in F}
        return frozenset(F)
    def L0(s): return [1 << (s.m + j) for j in range(s.m)]
    def level(s, F):
        L0 = s.frame(tuple(1 << j for j in range(2 * s.m)), s.L0())
        return s.m - (len(F & L0).bit_length() - 1)

def step(E, C, st):
    C = [dict(c) for c in C]
    t = st[0]
    def addto(d, k, v):
        w = d.get(k, 0) + v
        if abs(w) < 1e-15: d.pop(k, None)
        else: d[k] = w
    if t in ("move", "free"):
        r, M = st[1], st[2]
        C[r] = {(E.key(M @ E.mats[g]), l): c for (g, l), c in C[r].items()}
    elif t == "mix":
        r, summ = st[1], st[2]; new = {}
        for (sreg, mu, b) in summ:
            for (g, l), c in C[sreg].items():
                gk = g if b is None else E.key(b @ E.mats[g])
                addto(new, (gk, l), mu * c)
        C2 = [dict(c) for c in C]; C2[r] = new; C = C2
    elif t == "gate":
        R, A = st[1], st[2]; old = [dict(C[r]) for r in R]
        for i, r in enumerate(R):
            new = {}
            for j, sreg in enumerate(R):
                if A[i][j] != 0:
                    for k, c in old[j].items(): addto(new, k, A[i][j] * c)
            C[r] = new
    return C

def run(E, W, R, word):
    C = [dict() for _ in range(R)]
    for l in range(W): C[l] = {(E.IDk, l): 1.0}
    snaps = [C]
    for st in word: C = step(E, C, st); snaps.append(C)
    return snaps

def check_unit(E, W, snaps, deliver):
    """deliver: dict register -> input. Formal identity: content exactly u K y_l with |u| = 1."""
    C = snaps[-1]
    for r, l in deliver.items():
        if len(C[r]) != 1: return False
        (k, c), = C[r].items()
        if k != (E.Kk, l) or abs(abs(c) - 1) > 1e-12: return False
    return sorted(deliver.values()) == list(range(W))

def raw_arrays(E, W, R, word, deliver, f, rng):
    n = E.n; shape = (n,) * f
    def act(M, A):
        for ax in range(f): A = np.moveaxis(np.tensordot(M, A, axes=([1], [ax])), 0, ax)
        return A
    ys = [rng.standard_normal(shape) + 1j * rng.standard_normal(shape) for _ in range(W)]
    regs = [ys[l].copy() if l < W else np.zeros(shape, complex) for l in range(R)]
    for st in word:
        t = st[0]
        if t in ("move", "free"): regs[st[1]] = act(st[2], regs[st[1]])
        elif t == "mix":
            new = np.zeros(shape, complex)
            for (sreg, mu, b) in st[2]: new = new + mu * (regs[sreg] if b is None else act(b, regs[sreg]))
            regs[st[1]] = new
        elif t == "gate":
            old = [regs[r].copy() for r in st[1]]
            for i, r in enumerate(st[1]): regs[r] = sum(st[2][i][j] * old[j] for j in range(len(st[1])) if st[2][i][j] != 0)
    err = 0.0
    for r, l in deliver.items():
        target = act(E.Kmat, ys[l])
        # unit scalar u: fitted
        u = np.vdot(target, regs[r]) / np.vdot(target, target)
        err = max(err, float(np.max(np.abs(regs[r] - u * target))), abs(abs(u) - 1))
    return err

# ---------------- classes and liveness ----------------
def class_fn(E, mode, gauge=None):
    base = E.L0()
    if mode == "sp":      # Sp label g h_l
        return lambda g, l: E.comp(E.sp(g), gauge[l]) if gauge is not None else E.sp(g)
    if mode == "img":     # image frame gbar L0 (raw gauge)
        return lambda g, l: E.frame(E.sp(g), base)
    raise ValueError

def analyse(E, W, R, word, snaps, deliver, cls):
    T = len(word)
    nodes = [[set(cls(g, l) for (g, l) in snaps[t][r]) for r in range(R)] for t in range(T + 1)]
    live = [[set() for r in range(R)] for t in range(T + 1)]
    for r in deliver: live[T][r] = set(nodes[T][r])
    for t in range(T - 1, -1, -1):
        st = word[t]; ty = st[0]
        edges = []      # (r_from, kappa_from, r_to, kappa_to)
        touched = set()
        if ty in ("move", "free"):
            r = st[1]; touched = {r}
            for (g, l) in snaps[t][r]:
                gk = E.key(st[2] @ E.mats[g]); edges.append((r, cls(g, l), r, cls(gk, l)))
        elif ty == "mix":
            r = st[1]; touched = {r}
            for (sreg, mu, b) in st[2]:
                for (g, l) in snaps[t][sreg]:
                    gk = g if b is None else E.key(b @ E.mats[g])
                    k2 = cls(gk, l)
                    if k2 in nodes[t + 1][r]: edges.append((sreg, cls(g, l), r, k2))
        elif ty == "gate":
            touched = set(st[1])
            for i, r in enumerate(st[1]):
                for j, sreg in enumerate(st[1]):
                    if st[2][i][j] != 0:
                        for (g, l) in snaps[t][sreg]:
                            k2 = cls(g, l)
                            if k2 in nodes[t + 1][r]: edges.append((sreg, k2, r, k2))
        for r in range(R):
            if r not in touched:
                for k in nodes[t][r]:
                    if k in nodes[t + 1][r]: edges.append((r, k, r, k))
        for (r1, k1, r2, k2) in edges:
            if k2 in live[t + 1][r2]: live[t][r1].add(k1)
    return nodes, live

def moves_report(E, W, R, word, snaps, nodes, live, cls, levelf=None):
    """per paid move: (k_live, ascents of live classes, list of (level_before, level_after))."""
    out = []
    for t, st in enumerate(word):
        if st[0] != "move": continue
        r = st[1]; L = live[t][r]
        steps = []
        if levelf is not None:
            # map each live class to its image under the move (term-wise)
            seen = {}
            for (g, l) in snaps[t][r]:
                k = cls(g, l)
                if k in L and k not in seen: seen[k] = cls(E.key(st[2] @ E.mats[g]), l)
            steps = [(levelf(k), levelf(k2)) for k, k2 in seen.items()]
        out.append((len(L), steps))
    return out
