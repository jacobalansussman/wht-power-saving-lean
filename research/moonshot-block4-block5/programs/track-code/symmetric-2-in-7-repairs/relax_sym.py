"""relax_sym.py -- mw2-absorb4: the frame-sum relaxation matrix A(alpha) of a skeleton, with EXACT polynomial entries
(integer polynomials in the free content coefficients).  Phi([g] x_l) = e_(frame(g), l); A[:, i] = Phi(generator i).
Dumps A as a list of (row, generator, polynomial) and the row labels (frame distance, array).
Usage: relax_sym.py seed start [recipe.pkl]   (skeleton as in diag_border.py)"""
import sys, pickle
from relax_search import *
from collections import defaultdict

class Poly:
    """integer polynomial in nv variables: dict exponent-tuple -> int"""
    def __init__(self, d=None, nv=0):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}; self.nv = nv
    @staticmethod
    def const(c, nv):
        return Poly({(0,) * nv: c}, nv)
    @staticmethod
    def var(j, nv):
        e = [0] * nv; e[j] = 1
        return Poly({tuple(e): 1}, nv)
    def __add__(self, o):
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return Poly(d, self.nv)
    def __mul__(self, o):
        if not isinstance(o, Poly):
            return Poly({k: v * o for k, v in self.d.items()}, self.nv)
        d = defaultdict(int)
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                d[tuple(a + b for a, b in zip(k1, k2))] += v1 * v2
        return Poly(dict(d), self.nv)
    __rmul__ = __mul__
    def iszero(self):
        return not self.d
    def __repr__(self, names=None):
        if not self.d:
            return "0"
        out = []
        for k, v in sorted(self.d.items()):
            mon = "*".join((names[j] if names else f"a{j}") + (f"^{e}" if e > 1 else "") for j, e in enumerate(k) if e)
            out.append((str(v) if (v != 1 or not mon) else "") + ("*" if mon and v != 1 else "") + mon)
        return " + ".join(out)
    def ev(self, vals):
        return sum(v * np.prod([vals[j] ** e for j, e in enumerate(k)]) for k, v in self.d.items())

def sym_coefs(sk):
    """coefficient vectors of every generator as Poly in the free alphas (fixed alphas = 1)"""
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nv = len(free_idx)
    pos = {a: j for j, a in enumerate(free_idx)}
    al = [Poly.var(pos[a], nv) if a in pos else Poly.const(1, nv) for a in range(N)]
    out = []
    for i, (kind, *d) in enumerate(sk.coef_of):
        if kind == 'inp':
            out.append([Poly.const(1, nv)])
        elif kind == 'move':
            mv = sk.moves[d[0]]
            v = [Poly({}, nv) for _ in range(mv["n"])]
            for (a, src, si, ri, ch) in mv["maps"]:
                chr_ = complex(ch)
                assert abs(chr_.imag) < 1e-12 and abs(chr_.real - round(chr_.real)) < 1e-12
                for s_, r_ in zip(si, ri):
                    v[r_] = v[r_] + al[a] * out[src][s_] * int(round(chr_.real))
            out.append(v)
        else:
            out.append(out[sk.moves[d[0]]["rep"]])
    return out, free_idx

def relax_matrix(sk):
    G = sk.G
    cv, free_idx = sym_coefs(sk)
    fidx = {}; rowinfo = []
    for i, T in enumerate(sk.gens):
        for (l, g) in T:
            key = (G.frm(g), l)
            if key not in fidx:
                fidx[key] = len(fidx); rowinfo.append((G.dist(g), l, g))
    nv = len(free_idx)
    A = defaultdict(lambda: Poly({}, nv))
    for i, T in enumerate(sk.gens):
        for j, (l, g) in enumerate(T):
            r = fidx[(G.frm(g), l)]
            A[(r, i)] = A[(r, i)] + cv[i][j]
    A = {k: v for k, v in A.items() if not v.iszero()}
    Krow = fidx.get((G.frm(G.KID), 0))
    return A, rowinfo, Krow, free_idx

if __name__ == "__main__":
    m, W, H, shape, menu = 4, 2, 'Z2swap', 'P:inp,P:prev,L:prev,P:prev'.split(','), 2
    seed, s = int(sys.argv[1]), int(sys.argv[2])
    if len(sys.argv) > 3:
        d = pickle.load(open(sys.argv[3], 'rb'))
        sk = build_from_recipe(m, W, H, d['recipe'], d['skseed'])
    else:
        sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
    A, rowinfo, Krow, free_idx = relax_matrix(sk)
    names = [f"a{a}" for a in free_idx]
    print("free alphas", names, "rows", len(rowinfo), "K row", Krow, "generators", len(sk.gens))
    for r, (dd, l, g) in enumerate(rowinfo):
        ents = [(i, A[(r, i)]) for i in range(len(sk.gens)) if (r, i) in A]
        print(f"row {r:2d} dist {dd} arr {l}: " + " | ".join(f"g{i}: {p.__repr__(names)}" for i, p in ents))
    pickle.dump({"A": {k: v.d for k, v in A.items()}, "rowinfo": [(dd, l) for (dd, l, g) in rowinfo], "Krow": Krow,
                 "free_idx": free_idx, "ngen": len(sk.gens)}, open(f"relaxA_{seed}_{s}.pkl", "wb"))
