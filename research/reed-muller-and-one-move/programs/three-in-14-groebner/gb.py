"""gb.py: the exact algebra, theory-free and order-free.

For an instance and a set of forms S (for an excluded instance: the forms of its zero-block certificate), write the
polynomial system "the program's own equations at the forms in S" (ex2-families Lemma 1 / rm2-solve inst.py, with a
target per array):
    for every Q in S, with D = D_Q (moves where Q = 1) and P = arrays odd at Q:
        A_D N_D B_D^T = diag(P)        (W x W: a_(l,s) n_(s,t) b_(l',t) summed over s, t in D)
        N_D (I + G_D) = I              (N_D is the inverse of I + G_D)
in the unknowns a_(l,i), b_(l,i) (the columns of A and B at the moves i that occur), g_(i,j) for i != j (G has zero
diagonal and is otherwise FREE), n_(D; s,t).  A time-ordered program has G = H strictly upper triangular in its
order, and then N_D = (I + H_D)^-1 exists; so every program, in every time order, gives a solution of this system.
If 1 lies in the ideal (a Groebner basis over Q is {1}), no program exists in any time order, over C or any subfield.
For |D| = 1 the inverse is 1 and no n or g occurs.

Writes, per instance, Singular scripts (over Q with lift(I, 1) for the cofactors; over the prime 32003) and msolve
inputs (primes 65521 and 1073741827), all powers written x^k (there are none: every monomial is multilinear, which
gb.py asserts), and a JSON with the variable names and polynomials for verify_lift.py.
"""
import json, os, sys, itertools
from f2q import Field, Config, popc

PRIMES_MS = (65521, 1073741827)
PRIME_SING = 32003


class Sys:
    def __init__(self):
        self.vars = []; self.idx = {}; self.polys = []       # poly: dict monomial(tuple of var idx, sorted) -> int

    def v(self, name):
        if name not in self.idx:
            self.idx[name] = len(self.vars); self.vars.append(name)
        return self.idx[name]

    def add(self, poly):
        poly = {m: c for m, c in poly.items() if c}
        if poly: self.polys.append(poly)


def add_form(S, W, D, P, tag):
    """D: sorted list of move ids, P: set of odd arrays."""
    if len(D) == 1:
        i = D[0]
        for l in range(W):
            for l2 in range(W):
                poly = {tuple(sorted((S.v(f"a{l}_{i}"), S.v(f"b{l2}_{i}")))): 1}
                if l == l2 and l in P: poly[()] = poly.get((), 0) - 1
                S.add(poly)
        return
    n = {(s, t): S.v(f"n{tag}_{s}_{t}") for s in D for t in D}
    g = {(s, t): S.v(f"g{s}_{t}") for s in D for t in D if s != t}
    for s in D:                                   # N (I + G) = I
        for u in D:
            poly = {(n[(s, u)],): 1}
            for t in D:
                if t != u:
                    mono = tuple(sorted((n[(s, t)], g[(t, u)])))
                    poly[mono] = poly.get(mono, 0) + 1
            if s == u: poly[()] = -1
            S.add(poly)
    for l in range(W):                            # A N B^T = diag(P)
        for l2 in range(W):
            poly = {}
            for s in D:
                for t in D:
                    mono = tuple(sorted((S.v(f"a{l}_{s}"), n[(s, t)], S.v(f"b{l2}_{t}"))))
                    poly[mono] = poly.get(mono, 0) + 1
            if l == l2 and l in P: poly[()] = -1
            S.add(poly)


def build(cfg, forms):
    """forms: list of form bitmasks Q. Returns Sys."""
    S = Sys(); seenD = set()
    for Q in forms:
        Dm = cfg.Dmask[Q]; Pm = cfg.pmask[Q]
        if Dm in seenD: continue
        seenD.add(Dm)
        D = [i for i in range(cfg.M) if Dm >> i & 1]
        P = {l for l in range(cfg.W) if Pm >> l & 1}
        add_form(S, cfg.W, D, P, tag=f"D{Dm}")
    for p in S.polys:
        for mono in p: assert len(set(mono)) == len(mono), "non-multilinear monomial"
    return S


def poly_str(S, p, style):
    terms = []
    for mono, c in sorted(p.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        if style == "sing": vs = [f"x({i + 1})" for i in mono]
        else: vs = [f"x{i + 1}" for i in mono]
        if not vs: body = str(abs(c))
        elif abs(c) == 1: body = "*".join(vs)
        else: body = f"{abs(c)}*" + "*".join(vs)
        terms.append(("-" if c < 0 else "+") + body)
    s = "".join(terms)
    return s[1:] if s.startswith("+") else s


def write(S, base):
    NV = len(S.vars)
    polys_s = ",\n".join(poly_str(S, p, "sing") for p in S.polys)
    with open(base + "_Q.sing", "w") as f:
        f.write(f'ring r = 0, (x(1..{NV})), dp;\noption(redSB);\nshort = 0;\nideal I = {polys_s};\n'
                f'int t0 = timer;\nideal G = std(I);\n"GB-SIZE"; size(G);\n"GB-FIRST"; G[1];\n'
                f'if (G[1] == 1) {{\n  matrix T = lift(I, ideal(1));\n  "LIFT-BEGIN";\n  int i;\n'
                f'  for (i = 1; i <= size(I); i++) {{ string(T[i, 1]); }}\n  "LIFT-END";\n}}\n'
                f'"DIM"; dim(G);\nquit;\n')
    with open(base + "_p.sing", "w") as f:
        f.write(f'ring r = {PRIME_SING}, (x(1..{NV})), dp;\noption(redSB);\nideal I = {polys_s};\n'
                f'ideal G = std(I);\n"GB-SIZE"; size(G);\n"GB-FIRST"; G[1];\n"DIM"; dim(G);\nquit;\n')
    polys_m = ",\n".join(poly_str(S, p, "ms") for p in S.polys)
    for pr in PRIMES_MS:
        with open(base + f"_{pr}.ms", "w") as f:
            f.write(",".join(f"x{i + 1}" for i in range(NV)) + f"\n{pr}\n" + polys_m + "\n")
    json.dump(dict(vars=S.vars, polys=[[[list(m), c] for m, c in p.items()] for p in S.polys]),
              open(base + "_sys.json", "w"))
