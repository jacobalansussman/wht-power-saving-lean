# dualdfz.py -- the DUAL of DFZ (65) and its test (pf-duality, 14 Oct 2026).
# Linear polymatroid duality: for subspaces Y_1..Y_k (here the 7 letters), r*(S) = sum_{i in S} r(i) - r(E) + r(E\S)
# is again the rank function of subspaces over the same field (Z_i = (K cap W_{E\i})^perp in K*, K = kernel of
# (+)Y_i -> V; PROOF.md section 11).  Hence  sum_S c_S r*(S) >= 0  holds whenever sum_S c_S r(S) >= 0 does (DFZ (65):
# characteristic != 2).  The dual Fano matroid F7* (binary, rank 4) violates the dual inequality: it is genuinely
# 4-dimensional.  This file builds the dual term list (letters 0..6 = A,B,C,W,X,Y,Z) and tests it.
import random, itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

A_, B_, C_, W_, X_, Y_, Z_ = range(7)
def _Hs(c, *ls): return [(c, frozenset(ls)), (-c, frozenset())]
def _Hc(c, xs, ys): return [(c, frozenset(xs) | frozenset(ys)), (-c, frozenset(ys))]
_D = _Hs(-2, A_) + _Hs(-1, B_) + _Hs(-2, C_) + _Hs(1, W_) + _Hs(1, X_) + _Hs(1, Y_) + _Hs(1, Z_)
_D += _Hc(2, [A_], [Z_, Y_]) + _Hc(1, [B_], [X_, Z_]) + _Hc(2, [C_], [A_, X_]) + _Hc(3, [X_], [W_, Y_]) + _Hc(3, [Z_], [W_, C_])
_D += _Hc(5, [W_], [A_, B_]) + _Hc(5, [Y_], [B_, C_]) + _Hs(5, A_) + _Hs(5, B_) + _Hs(5, C_) + _Hs(-5, A_, B_, C_)
E7 = frozenset(range(7))
def dualize(terms):
    acc = {}
    for c, S in terms:
        # r*(S) = sum_{i in S} r({i}) - r(E) + r(E \ S)
        for i in S: acc[frozenset({i})] = acc.get(frozenset({i}), 0) + c
        acc[E7] = acc.get(E7, 0) - c
        acc[E7 - S] = acc.get(E7 - S, 0) + c
    out = [(c, s) for s, c in acc.items() if c != 0 and s]
    tot = sum(c for c, s in out)
    if tot: out.append((-tot, frozenset()))   # r_T(S) = r(T u S) - r(T): the -r(T) parts, needed when T is non-empty
    return out
DFZD = dualize(_D)          # each term (c, S) means c * r(T u S); sum >= 0 for subspaces over fields with 2 invertible

def rank_mod(rows, p):
    rows = [r[:] for r in rows if any(x % p for x in r)]; rk = 0
    if not rows: return 0
    for c in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p - 2, p) if p > 2 else 1
        rows[rk] = [(v * inv) % p for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                f = rows[i][c]; rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def value(terms, sp, p):
    return sum(c * rank_mod([v for l in s for v in sp[l]], p) for c, s in terms)

if __name__ == "__main__":
    # F7 (DFZ's labelling) and its Gale dual F7*: B = basis of the orthogonal complement of the row space of A over GF(2)
    cols = [1, 2, 4, 3, 5, 6, 7]                               # A,B,C,W,X,Y,Z as vectors in F_2^3
    Arows = [[(c >> b) & 1 for c in cols] for b in range(3)]   # 3 x 7
    null = [x for x in itertools.product(range(2), repeat=7) if all(sum(a * b for a, b in zip(r, x)) % 2 == 0 for r in Arows)]
    basis = []
    for x in null:
        if not any(x):
            continue
        span = {tuple([0] * 7)}
        for b in basis:
            span |= {tuple((u + v) % 2 for u, v in zip(s, b)) for s in span}
        if tuple(x) not in span: basis.append(list(x))
    assert len(basis) == 4
    dual_vecs = [[basis[r][j] for r in range(4)] for j in range(7)]  # column j of the 4 x 7 matrix
    f7 = {l: [[(cols[l] >> b) & 1 for b in range(3)]] for l in range(7)}
    f7s = {l: [dual_vecs[l]] for l in range(7)}
    print("F7  over GF(2): DFZ  =", value(_D, f7, 2), " dual-DFZ =", value(DFZD, f7, 2))
    print("F7* over GF(2): DFZ  =", value(_D, f7s, 2), " dual-DFZ =", value(DFZD, f7s, 2), " (dual-DFZ must be < 0)")
    print("F7* vectors (A..Z):", [''.join(map(str, v)) for v in dual_vecs])
    rng = random.Random(3); worst = None; n = 0
    for p in (3, 5, 7, 1000003):
        for t in range(3000 if p < 100 else 1000):
            dim = rng.randint(2, 6); mode = t % 3
            base = [[rng.randrange(p) for _ in range(dim)] for _ in range(4)]
            sp = {}
            for l in range(7):
                k = rng.randint(0, 2)
                vecs = [[rng.randrange(p) for _ in range(dim)] for _ in range(k)]
                if mode:
                    vecs = [[sum(rng.choice([0, 1, 1, p - 1]) * b[i] for b in base) % p for i in range(dim)] for _ in range(max(k, 1))]
                sp[l] = vecs if vecs else [[0] * dim]
            v = value(DFZD, sp, p); n += 1
            worst = v if worst is None else min(worst, v)
            assert v >= 0, ("VIOLATION", p, sp, v)
        # also the F7* vectors read over GF(p): must satisfy it
        assert value(DFZD, {l: [dual_vecs[l]] for l in range(7)}, p) >= 0
    print("dual-DFZ on", n, "random arrangements over GF(3), GF(5), GF(7), GF(1000003): min", worst, "(no violation)")
