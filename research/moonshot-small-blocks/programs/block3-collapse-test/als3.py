"""als3.py -- ms-xw3.  Block-coordinate exact least squares (ALS) over a DAG of paid moves at block 3.

A program: moves i = 0..t-1, y_i = C_i D_{z_i} C_i^-1 ( content_i ), outputs out_l, l < W.
A content / output is a list of TERMS:
   ('in', word, l, c)          fixed: c * word(v_l)
   ('y',  word, j, c)          fixed: c * word(y_j)
   ('Gin', key, menu)          GENERAL: sum over b in menu, l < W of u[b,l] * b(v_l)
   ('Gy',  key, menu, j)       GENERAL: sum over b in menu of u[b] * b(y_j)
Every GENERAL slot has a key; slots with the same key[0] form one block.  BCD: for each block in turn (all other
blocks fixed at their current numbers) the outputs are affine in the block's unknowns provided no path crosses two
slots of the block (asserted); that linear least squares is solved exactly, so the residual never increases.
Values are carried as (constant part, linear part or None), so only quantities that depend on the active block
carry unknown columns.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, flat


class Val:
    __slots__ = ("c", "L")

    def __init__(self, c, L=None):
        self.c = c; self.L = L          # c: f -> (n, S); L: f -> (n, S, U) or None


def vmap(G, fn, X):
    return Val(fn(X.c), None if X.L is None else fn(X.L))


def vadd(X, Y, coef=1.0):
    c = {f: X.c[f] + coef * Y.c[f] for f in X.c}
    if X.L is None and Y.L is None: L = None
    elif X.L is None: L = {f: coef * Y.L[f] for f in Y.L}
    elif Y.L is None: L = X.L
    else: L = {f: X.L[f] + coef * Y.L[f] for f in X.L}
    return Val(c, L)


class Prog:
    def __init__(self, G, moves, outs, W=2):
        self.G = G; self.moves = moves; self.outs = outs; self.W = W
        self.slots = {}
        for mv in moves:
            for t in mv["content"]:
                if t[0] in ("Gin", "Gy"): self._reg(t)
        for ol in outs:
            for t in ol:
                if t[0] in ("Gin", "Gy"): self._reg(t)
        self.blocks = sorted(set(k[0] for k in self.slots))
        self.img_in = {}

    def _reg(self, t):
        key = t[1]
        nunk = len(t[2]) * (self.W if t[0] == "Gin" else 1)
        self.slots[key] = (t[0], t[2], nunk)

    def inimg(self, b, l):
        if (b, l) not in self.img_in:
            self.img_in[(b, l)] = self.G.W_(b, self.G.inp(l))
        return self.img_in[(b, l)]

    def sparse_op(self, menu, u, f):
        import scipy.sparse as sp
        G = self.G; n = G.n[f]
        rows = []; cols = []; data = []
        PHc = np.array([1, 1j, -1, -1j])
        for bi, b in enumerate(menu):
            if u[bi] == 0: continue
            g = G.word(b)[f]
            rows.append(g.sig); cols.append(np.arange(n)); data.append(u[bi] * PHc[g.ph])
        if not rows: return sp.csr_matrix((n, n), dtype=complex)
        return sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))

    def evaluate(self, vals, active=None):
        G = self.G; FS = G.FS
        keys = [k for k in sorted(self.slots, key=str) if k[0] == active] if active is not None else []
        U = sum(self.slots[k][2] for k in keys)
        start = {}; pos = 0
        for k in keys: start[k] = pos; pos += self.slots[k][2]
        zc = lambda: {f: np.zeros((G.n[f], G.S), complex) for f in FS}

        def term_val(t, ys):
            if t[0] == "in":
                _, w, l, c = t; return Val({f: c * self.inimg(w, l)[f] for f in FS})
            if t[0] == "y":
                _, w, j, c = t; X = vmap(G, lambda A: G.W_(w, A), ys[j])
                return Val({f: c * X.c[f] for f in FS}, None if X.L is None else {f: c * X.L[f] for f in FS})
            key, menu = t[1], t[2]
            if key in start:
                L = {f: np.zeros((G.n[f], G.S, U), complex) for f in FS}
                k0 = start[key]
                if t[0] == "Gin":
                    for bi, b in enumerate(menu):
                        for l in range(self.W):
                            im = self.inimg(b, l)
                            for f in FS: L[f][:, :, k0 + bi * self.W + l] = im[f]
                else:
                    j = t[3]; assert ys[j].L is None, "bilinear: two active slots on one path"
                    for bi, b in enumerate(menu):
                        im = G.W_(b, ys[j].c)
                        for f in FS: L[f][:, :, k0 + bi] = im[f]
                return Val(zc(), L)
            u = vals[key]; acc = Val(zc())
            if t[0] == "Gin":
                for bi, b in enumerate(menu):
                    for l in range(self.W):
                        c = u[bi * self.W + l]
                        if c != 0: acc = vadd(acc, Val(self.inimg(b, l)), c)
            else:
                j = t[3]
                if ys[j].L is None:
                    for bi, b in enumerate(menu):
                        if u[bi] != 0: acc = vadd(acc, Val(G.W_(b, ys[j].c)), u[bi])
                else:                      # sparse operator sum_b u_b b applied to constant and linear parts at once
                    c = {}; L = {}
                    for f in FS:
                        Sm = self.sparse_op(menu, u, f)
                        c[f] = Sm @ ys[j].c[f]
                        n, S_, U_ = ys[j].L[f].shape
                        L[f] = (Sm @ ys[j].L[f].reshape(n, -1)).reshape(n, S_, U_)
                    acc = Val(c, L)
            return acc

        ys = []
        for mv in self.moves:
            cont = Val(zc())
            for t in mv["content"]: cont = vadd(cont, term_val(t, ys))
            ys.append(vmap(G, lambda A, mv=mv: G.move(mv["z"], mv.get("C", "1"), A), cont))
        outs = []
        for ol in self.outs:
            o = Val(zc())
            for t in ol: o = vadd(o, term_val(t, ys))
            outs.append(o)
        return outs, U, ys

    def init_vals(self, rng, scale=0.3, density=1.0):
        v = {}
        for key, (kind, menu, nunk) in self.slots.items():
            x = scale * (rng.standard_normal(nunk) + 1j * rng.standard_normal(nunk)) / np.sqrt(2)
            if density < 1.0: x = x * (rng.random(nunk) < density)
            v[key] = x
        return v


def targets_K(G, W):
    return [G.K(G.inp(l)) for l in range(W)]


def fl(X):
    return np.concatenate([X[f].reshape(-1) for f in sorted(X)])


def residual(P, vals, T):
    outs, U, _ = P.evaluate(vals, None)
    num = sum(np.linalg.norm(fl(o.c) - fl(T[l])) ** 2 for l, o in enumerate(outs))
    den = sum(np.linalg.norm(fl(T[l])) ** 2 for l in range(len(T)))
    return float(np.sqrt(num / den))


def per_output(P, vals, T):
    outs, U, _ = P.evaluate(vals, None)
    return [float(np.linalg.norm(fl(o.c) - fl(T[l])) / np.linalg.norm(fl(T[l]))) for l, o in enumerate(outs)]


def bcd(P, T, vals, sweeps=30, log=None, tol=1e-10, stall=0.995, order=None, giveup=None):
    den = sum(np.linalg.norm(fl(T[l])) ** 2 for l in range(len(T)))
    hist = []
    blocks = order or P.blocks
    res = None
    for sw in range(sweeps):
        for blk in blocks:
            if blk == "R":                      # readouts: block diagonal over outputs -> solve each output alone
                outs, U, _ = P.evaluate(vals, blk)
                keys = [k for k in sorted(P.slots, key=str) if k[0] == "R"]
                offs = {}; pos = 0
                for k in keys: offs[k] = (pos, pos + P.slots[k][2]); pos += P.slots[k][2]
                tot = 0.0
                for l, o in enumerate(outs):
                    kl = [k for k in keys if k[1] == l]
                    idx = np.concatenate([np.arange(*offs[k]) for k in kl])
                    A = np.concatenate([o.L[f].reshape(-1, U)[:, idx] for f in sorted(o.L)], axis=0)
                    b = fl(T[l]) - fl(o.c)
                    x = _lstsq(A, b)
                    tot += np.linalg.norm(A @ x - b) ** 2
                    for k in kl:
                        n = P.slots[k][2]; vals[k] = x[:n]; x = x[n:]
                res = float(np.sqrt(tot / den))
                continue
            outs, U, _ = P.evaluate(vals, blk)
            if U == 0: continue
            A = np.concatenate([np.concatenate([o.L[f].reshape(-1, U) for f in sorted(o.L)], axis=0)
                                if o.L is not None else np.zeros((fl(o.c).size, U), complex) for o in outs], axis=0)
            b = np.concatenate([fl(T[l]) - fl(o.c) for l, o in enumerate(outs)])
            x = _lstsq(A, b)
            res = float(np.sqrt(np.linalg.norm(A @ x - b) ** 2 / den))
            for key in sorted(P.slots, key=str):
                if key[0] == blk:
                    n = P.slots[key][2]
                    vals[key] = x[:n]; x = x[n:]
        hist.append(res)
        if log: log(f"    sweep {sw}: {res:.4e}")
        if res < tol: break
        if sw >= 6 and hist[-1] > stall * hist[-5]: break
        if giveup is not None and sw + 1 >= giveup[0] and hist[-1] > giveup[1]: break
    return hist[-1], hist


def _lstsq(A, b):
    try:
        return np.linalg.lstsq(A, b, rcond=None)[0]
    except np.linalg.LinAlgError:              # rare gelsd failure: fall back to QR with pivoting
        import scipy.linalg as sl
        return sl.lstsq(A, b, lapack_driver="gelsy", check_finite=False)[0]
