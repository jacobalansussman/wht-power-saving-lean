#!/usr/bin/env python3
"""vgame.py (v-pk2-orbits, skeptic): an independent frame game C[Sp(2m,2)] and exact closure solver.
Written from the definitions, not from fpa.py.  Group elements are tuples of COLUMN images (g[j] = g(e_j)), bits
0..m-1 = x-part a, bits m..2m-1 = z-part b.  omega((a|b),(a'|b')) = a.b' + b.a'.  T_u(x) = x + omega(x,u) u.
K = prod_i T_(e_i|0) : (a|b) -> (a+b|b).  L0 = {(0|b)}.  frame(g) = g^-1 L0.  level(g) = m - dim(L0 cap g^-1 L0).
P = {g : g L0 = L0}.  Contents: dict (input l, element g) -> Fraction.  Free steps multiply on the left by C[P]."""
from fractions import Fraction as Fr

def pc(x): return bin(x).count("1") & 1

class G:
    def __init__(s, m):
        s.m = m; s.n = 2 * m; s.ma = (1 << m) - 1
        s.I = tuple(1 << j for j in range(s.n))
        s.L0 = [1 << (m + j) for j in range(m)]
        s.K = s.I
        for i in range(m): s.K = s.comp(s.T(1 << i), s.K)
        s._fr = {}
    def ap(s, g, v):
        out = 0; j = 0
        while v:
            if v & 1: out ^= g[j]
            v >>= 1; j += 1
        return out
    def comp(s, g, h): return tuple(s.ap(g, h[j]) for j in range(s.n))
    def om(s, v, w):
        m = s.m
        return pc(((v & s.ma) & (w >> m)) ^ ((v >> m) & (w & s.ma)))
    def T(s, u): return tuple((1 << j) ^ (u if s.om(1 << j, u) else 0) for j in range(s.n))
    def inv(s, g):
        # generic GF(2) inversion: solve g x = e_j by Gaussian elimination on [g | I] (rows = equations)
        n = s.n
        rows = []
        for i in range(n):            # row i of matrix g: bit j = (g[j] >> i) & 1
            r = 0
            for j in range(n):
                if (g[j] >> i) & 1: r |= 1 << j
            rows.append((r, 1 << i))
        for c in range(n):
            p = next(k for k in range(c, n) if (rows[k][0] >> c) & 1)
            rows[c], rows[p] = rows[p], rows[c]
            for k in range(n):
                if k != c and (rows[k][0] >> c) & 1:
                    rows[k] = (rows[k][0] ^ rows[c][0], rows[k][1] ^ rows[c][1])
        # now rows[i] = (e_i, row i of g^-1)
        cols = []
        for j in range(n):
            v = 0
            for i in range(n):
                if (rows[i][1] >> j) & 1: v |= 1 << i
            cols.append(v)
        return tuple(cols)
    def is_symp(s, g):
        return all(s.om(g[i], g[j]) == s.om(1 << i, 1 << j) for i in range(s.n) for j in range(s.n))
    @staticmethod
    def span_rref(vecs):
        basis = []
        for v in vecs:
            for b in basis:
                if v ^ b < v: v ^= b
            if v:
                basis = [b ^ v if b ^ v < b else b for b in basis]
                basis.append(v)
        # canonical reduced form
        basis.sort(reverse=True)
        red = []
        for v in basis:
            for b in red:
                if v ^ b < v: v ^= b
            red.append(v)
        out = []
        red.sort(reverse=True)
        for i, v in enumerate(red):
            for j in range(len(red)):
                if j != i:
                    hb = red[j].bit_length() - 1
                    if (v >> hb) & 1: v ^= red[j]
            out.append(v)
        return tuple(sorted(out))
    def dimsp(s, vecs): return len(s.span_rref(vecs))
    def frame(s, g):
        if g not in s._fr:
            gi = s.inv(g)
            s._fr[g] = s.span_rref([s.ap(gi, v) for v in s.L0])
        return s._fr[g]
    def dist(s, X, Y):  # Lagrangians given as bases
        return s.m - (2 * s.m - s.dimsp(list(X) + list(Y)))
    def level(s, g): return s.dist(s.span_rref(s.L0), s.frame(g))
    def inP(s, g): return s.span_rref([s.ap(g, v) for v in s.L0]) == s.span_rref(s.L0)

def lm(Gp, h, c):
    out = {}
    for (l, g), v in c.items():
        k = (l, Gp.comp(h, g)); out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}
def add(*pairs):
    out = {}
    for a, c in pairs:
        for k, v in c.items(): out[k] = out.get(k, 0) + a * v
    return {k: v for k, v in out.items() if v != 0}

def gauss_solve(eqs, rhs):
    """eqs: list of dict col -> Fraction; rhs list. Returns (consistent, solution dict)."""
    piv = {}; order = []
    for e, b in zip(eqs, rhs):
        e = dict(e); b = Fr(b)
        changed = True
        while True:
            hit = [c for c in e if c in piv]
            if not hit: break
            c = hit[0]; f = e[c]; pe, pb = piv[c]
            for k, v in pe.items():
                nv = e.get(k, 0) - f * v
                if nv == 0: e.pop(k, None)
                else: e[k] = nv
            b -= f * pb
        if not e:
            if b != 0: return False, None
            continue
        c = max(e); f = e[c]
        e = {k: v / f for k, v in e.items()}; b /= f
        piv[c] = (e, b); order.append(c)
    # back substitution: free vars 0; resolve in reverse order of creation
    sol = {}
    for c in reversed(order):
        pe, pb = piv[c]
        val = pb - sum(v * sol.get(k, 0) for k, v in pe.items() if k != c)
        sol[c] = val
    return True, sol

def closure_solve(Gp, states, target, cap=30000):
    """exact: is target in the C[P]-span of the states?  Closure over (state, p) columns from the target rows."""
    byframe = {}
    for r, st in enumerate(states):
        for (l, g), v in st.items():
            byframe.setdefault((l, Gp.frame(g)), []).append((r, g))
    invc = {}
    def inv(g):
        if g not in invc: invc[g] = Gp.inv(g)
        return invc[g]
    rows = set(target); queue = list(target); cols = {}
    while queue:
        if len(cols) > cap: return "cap", None
        (l, g) = queue.pop()
        for (r, h) in byframe.get((l, Gp.frame(g)), []):
            p = Gp.comp(g, inv(h))
            if (r, p) in cols: continue
            col = lm(Gp, p, states[r])
            cols[(r, p)] = col
            for k in col:
                if k not in rows: rows.add(k); queue.append(k)
    ckeys = list(cols)
    eqd = {}
    for ci, ck in enumerate(ckeys):
        for k, v in cols[ck].items(): eqd.setdefault(k, {})[ci] = v
    rl = list(rows)
    ok, sol = gauss_solve([eqd.get(k, {}) for k in rl], [target.get(k, 0) for k in rl])
    if not ok: return "no", None
    # verify
    acc = {}
    for ci, x in sol.items():
        if x: acc = add((1, acc), (x, cols[ckeys[ci]]))
    assert add((1, acc), (-1, target)) == {}, "solution failed verification"
    return "yes", [(ckeys[ci], x) for ci, x in sol.items() if x]

def psi_screen(Gp, states, l):
    """[L1] e_l in the scalar span of the frame profiles of the states"""
    tgt = (l, Gp.frame(Gp.K))
    vecs = []
    for st in states:
        d = {}
        for (ll, g), v in st.items():
            k = (ll, Gp.frame(g)); d[k] = d.get(k, 0) + v
        vecs.append({k: v for k, v in d.items() if v != 0})
    rows = set(k for v in vecs for k in v) | {tgt}
    eqd = {}
    for ci, v in enumerate(vecs):
        for k, x in v.items(): eqd.setdefault(k, {})[ci] = x
    rl = list(rows)
    ok, _ = gauss_solve([eqd.get(k, {}) for k in rl], [1 if k == tgt else 0 for k in rl])
    return ok
