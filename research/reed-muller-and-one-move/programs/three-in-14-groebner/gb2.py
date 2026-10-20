"""gb2.py: the exact algebra, theory-free, one time order at a time (the program's own equations, as in
rm2-solve's inst.py / tri_direct.py, restricted to a set U of moves).

For an instance, a move set U and a time order pi of U:
  unknowns a_(l,i), b_(l,i) (l < W, i in U) and h_(i,j) for i before j in pi;
  for EVERY form Q whose support D_Q lies inside U (D_Q nonempty), with P = arrays odd at Q:
     F_Q[l, l'] = sum over chains s_1 < ... < s_k in D_Q (in the order pi) of
                  (-1)^(k-1) a_(l,s_1) h_(s_1,s_2) ... h_(s_(k-1),s_k) b_(l',s_k)  =  [l = l' and l in P].
These are a subset of the full system's equations (every form of the full system with support in U), and they
involve only the unknowns at moves of U.  So if 1 is in the ideal for EVERY relative order of U, the instance has
no program in any time order (over C, hence over any subfield).  Orders that differ only by swapping two copies of
the same direction give the same system up to renaming, so one order per multiset permutation of the directions; and
reversing the order while swapping A with B and transposing H maps F_Q to its transpose, and the targets are
diagonal, so an order and its reverse give the same answer: one of each pair is kept.

All monomials are multilinear (asserted); polynomials are written with explicit '*' and no powers.
"""
import itertools, json
from f2q import popc


def supports_in(cfg, U):
    """dict support (bitmask inside U) -> parity bitmask, for every form with nonempty support inside U."""
    Um = sum(1 << i for i in U); out = {}
    for Q in range(1, len(cfg.Dmask)):
        d = cfg.Dmask[Q]
        if d and (d & ~Um) == 0:
            p = cfg.pmask[Q]
            if d in out: assert out[d] == p, "span violated"
            out[d] = p
    return out


def orders(cfg, U):
    """one representative order per multiset permutation of the directions of U."""
    dirs = [cfg.moves[i] for i in U]
    copies = {}
    for i in U: copies.setdefault(cfg.moves[i], []).append(i)
    seen = set()
    for perm in itertools.permutations(range(len(U))):
        seq = tuple(dirs[p] for p in perm)
        if seq in seen: continue
        seen.add(seq)
        if seq[::-1] < seq: continue      # time reversal (A <-> B, H -> H^T) maps order seq to its reverse
        used = {z: 0 for z in copies}; order = []
        for z in seq:
            order.append(copies[z][used[z]]); used[z] += 1
        yield order


class Ring:
    """variable names shared by all orders of one instance: a, b at moves of U; h for every ordered pair."""
    def __init__(self, W, U):
        self.names = [f"a{l}_{i}" for l in range(W) for i in U] + [f"b{l}_{i}" for l in range(W) for i in U] + \
                     [f"h{i}_{j}" for i in U for j in U if i != j]
        self.idx = {n: k for k, n in enumerate(self.names)}


def system(cfg, R, order, sup):
    """list of polys (dict mono(tuple of var idx) -> int) for one time order."""
    pos = {i: k for k, i in enumerate(order)}; W = cfg.W
    polys = []
    for d, p in sorted(sup.items()):
        S = sorted((i for i in order if d >> i & 1), key=lambda i: pos[i])
        chains = []
        for k in range(1, len(S) + 1):
            for ch in itertools.combinations(S, k):        # combinations keep the time order of S
                chains.append(ch)
        for l in range(W):
            for l2 in range(W):
                poly = {}
                for ch in chains:
                    vs = [R.idx[f"a{l}_{ch[0]}"]] + [R.idx[f"h{ch[t]}_{ch[t + 1]}"] for t in range(len(ch) - 1)] + \
                         [R.idx[f"b{l2}_{ch[-1]}"]]
                    mono = tuple(sorted(vs)); assert len(set(mono)) == len(mono)
                    poly[mono] = poly.get(mono, 0) + (-1) ** (len(ch) - 1)
                if l == l2 and (p >> l & 1): poly[()] = poly.get((), 0) - 1
                poly = {m: c for m, c in poly.items() if c}
                if poly: polys.append(poly)
    return polys


def pstr(p, style):
    terms = []
    for mono, c in sorted(p.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        vs = [f"x({i + 1})" if style == "sing" else f"x{i + 1}" for i in mono]
        body = str(abs(c)) if not vs else ("*".join(vs) if abs(c) == 1 else f"{abs(c)}*" + "*".join(vs))
        terms.append(("-" if c < 0 else "+") + body)
    s = "".join(terms)
    return s[1:] if s.startswith("+") else s


def used_vars(polys):
    return sorted({v for p in polys for m in p for v in m})


def msolve_text(polys, prime):
    """msolve input with only the variables that occur (renumbered)."""
    uv = used_vars(polys); ren = {v: k for k, v in enumerate(uv)}
    rp = [{tuple(sorted(ren[v] for v in m)): c for m, c in p.items()} for p in polys]
    return ",".join(f"x{k + 1}" for k in range(len(uv))) + f"\n{prime}\n" + ",\n".join(pstr(p, "ms") for p in rp) + "\n"


def singular_text(R, systems, char, label, lift=False):
    """one Singular script: one ring, one ideal per order; prints 'ORDER k: UNIT' or 'ORDER k: NOT-UNIT dim d'.
    With lift=True, for every UNIT order also the cofactors of lift(I, 1) between 'LIFT k BEGIN' and 'LIFT k END'."""
    NV = len(R.names)
    out = [f'ring r = {char}, (x(1..{NV})), dp;', 'option(redSB);', 'short = 0;', 'int nunit = 0; int nnot = 0;',
           'ideal I; ideal G; matrix T; int i;']
    for k, polys in enumerate(systems):
        out.append("I = " + ",\n".join(pstr(p, "sing") for p in polys) + ";")
        lt = (f' T = lift(I, ideal(1)); "LIFT {k} BEGIN"; for (i = 1; i <= ncols(I); i++) {{ string(T[i, 1]); }} "LIFT {k} END";'
              if lift else "")
        out.append(f'G = std(I); if (G[1] == 1) {{ nunit++; "ORDER {k}: UNIT";{lt} }} else {{ nnot++; "ORDER {k}: NOT-UNIT dim " + string(dim(G)); }}')
    out.append(f'"SUMMARY {label} char {char}: orders " + string({len(systems)}) + " unit " + string(nunit) + " not-unit " + string(nnot);')
    out.append("quit;")
    return "\n".join(out) + "\n"
