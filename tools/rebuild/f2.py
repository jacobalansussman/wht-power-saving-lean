"""f2.py: subspaces of F_2^h as tuples of integer masks.  Library of the rebuild scripts (tools/rebuild/).

Canonical form of a subspace = its reduced row echelon basis, pivots = leading bits, sorted by pivot descending
(the reduced echelon form is unique, so equal subspaces have equal tuples).
Form on F_2^h: the standard dot product a.b = parity(a & b).  "unit vector" = odd weight (a.a = 1).
Self-test: python3 -I -B f2.py
Two functions carry ADAPTED code: perp() follows perp() statement by statement, and the second half of rref()
(from "for q, y in list(piv.items())" to the return line) follows the second half of basis() with other names;
the two loops of rref() that reduce a row are written differently.  Both outside functions are in
scripts/paired_cube/frames.py of the outside repository (github.com/CrocSwap/integer-mult-bounds, Apache-2.0,
commit 4a3c769; read as text, never run; "Copyright 2026 icekylinx. Apache-2.0.").  NOTICE, section 6.
"""


def rref(rows):
    piv = {}
    for x in rows:
        while x:
            p = x.bit_length() - 1
            y = piv.get(p)
            if y is None:
                break
            x ^= y
        if x:
            p = x.bit_length() - 1
            # reduce: clear lower pivots from x, then clear bit p from the others
            for q, y in piv.items():
                if q < p and (x >> q) & 1:
                    x ^= y
            for q, y in list(piv.items()):
                if (y >> p) & 1:
                    piv[q] = y ^ x
            piv[p] = x
    return tuple(piv[p] for p in sorted(piv, reverse=True))


def perp(rows, h):
    """orthogonal complement (dot product) of span(rows) in F_2^h, canonical"""
    rows = rref(rows)
    piv = {r.bit_length() - 1: r for r in rows}
    out = []
    for j in range(h):
        if j in piv:
            continue
        x = 1 << j
        for p, r in piv.items():
            if (r >> j) & 1:
                x |= 1 << p
        out.append(x)
    return rref(out)


def inside(A, B):
    """span(A) <= span(B), B canonical (rref)"""
    pv = [(y.bit_length() - 1, y) for y in B]
    for x in A:
        for p, y in pv:
            if (x >> p) & 1:
                x ^= y
        if x:
            return False
    return True


def dot(a, b):
    return bin(a & b).count("1") & 1


def radical(U, h):
    """U cap U^perp, canonical"""
    U = list(U)
    if not U:
        return ()
    # kernel of the Gram matrix: combinations c with sum c_i U_i orthogonal to all U_j
    k = len(U)
    rows = []
    for i in range(k):
        g = 0
        for j in range(k):
            if dot(U[i], U[j]):
                g |= 1 << j
        rows.append((g, 1 << i))
    # eliminate on g, track combination
    piv = {}
    ker = []
    for g, c in rows:
        while g:
            p = g.bit_length() - 1
            if p in piv:
                g2, c2 = piv[p]
                g ^= g2
                c ^= c2
            else:
                piv[p] = (g, c)
                break
        if not g:
            ker.append(c)
    out = []
    for c in ker:
        x = 0
        for i in range(k):
            if (c >> i) & 1:
                x ^= U[i]
        out.append(x)
    return rref(out)


def kind(U, h):
    """classification of a subspace as a frame of the label calculus of this repository:
    'zero' / 'label' (nondegenerate and contains a unit vector: has an orthonormal basis) /
    'alt' (nondegenerate, every vector even: no unit vector) / 'deg' (U cap U^perp != 0)."""
    if not U:
        return "zero"
    if radical(U, h):
        return "deg"
    if all(bin(x).count("1") % 2 == 0 for x in U):
        return "alt"
    return "label"


def step_kind(A, B, h):
    """a block = one role going from frame A up to frame B (A <= B, A != B).
    'legal'  : A and B are labels (or zero) and the residual B cap A^perp is a label, so the step is an
               orthonormal family of unit vectors (rule L1 of the carrier specification, see tools/ccheck.py);
    'alt'    : A, B nondegenerate but the residual (or an end) is alternating: no orthonormal family;
    'deg'    : A or B is degenerate."""
    ka, kb = kind(A, h), kind(B, h)
    if ka == "deg" or kb == "deg":
        return "deg"
    if ka == "alt" or kb == "alt":
        return "alt"
    # residual = B cap A^perp ; A nondegenerate so B = A (+) residual
    if A:
        Bp = perp(B, h)
        res = perp(tuple(Bp) + tuple(A), h)   # (B^perp + A)^perp = B cap A^perp
    else:
        res = B
    kr = kind(res, h)
    assert kr != "deg" and len(res) == len(B) - len(A)
    return "legal" if kr == "label" else "alt"


def orthonormal(Rb):
    """orthonormal basis (unit vectors, pairwise orthogonal) of span(Rb), which must be nondegenerate and contain
    a unit vector.  Over F_2 the map x -> x.x is linear (parity of the weight)."""
    vecs, out = list(Rb), []
    while vecs:
        j = next((i for i, w in enumerate(vecs) if bin(w).count("1") & 1), None)
        if j is not None:
            u = vecs.pop(j)
            vecs = [w ^ u if dot(w, u) else w for w in vecs]
            out.append(u)
        else:                       # the rest is alternating: trade the last unit vector with a hyperbolic pair
            assert out, "alternating: no orthonormal basis"
            e, a = out.pop(), vecs.pop(0)
            j = next(i for i, w in enumerate(vecs) if dot(a, w))
            b = vecs.pop(j)
            vecs = [w ^ (a if dot(w, b) else 0) ^ (b if dot(w, a) else 0) for w in vecs]
            out += [e ^ a, e ^ b, e ^ a ^ b]
    assert all(dot(x, y) == (1 if i == j else 0) for i, x in enumerate(out) for j, y in enumerate(out))
    assert rref(out) == rref(Rb)
    return out


def step_dirs(A, B, h):
    """the new directions of the block A -> B: (kind, dirs) with span(A + dirs) = B and len(dirs) = dim B - dim A.
    legal: an ORTHONORMAL family orthogonal to A (exactly a carrier-cert/1 `dirs` list);
    alt / deg: any independent completion (here: the vectors of B's normal-form basis that are needed)."""
    k = step_kind(A, B, h)
    if k == "legal":
        res = perp(tuple(perp(B, h)) + tuple(A), h) if A else B
        return k, orthonormal(res)
    cur, dirs = tuple(A), []
    for b in B:
        if not inside((b,), cur):
            dirs.append(b)
            cur = rref(cur + (b,))
    assert cur == tuple(B) and len(dirs) == len(B) - len(A)
    return k, dirs


def selftest():
    import random
    rnd = random.Random(7)
    h = 10
    for _ in range(300):
        rows = [rnd.randrange(1, 1 << h) for _ in range(rnd.randrange(1, 7))]
        U = rref(rows)
        P = perp(U, h)
        assert len(U) + len(P) == h and all(dot(a, b) == 0 for a in U for b in P)
        assert perp(P, h) == U and inside(rows, U)
        rad = radical(U, h)
        assert inside(rad, U) and inside(rad, P)
        # brute force radical dimension
        cnt = 0
        for c in range(1 << len(U)):
            x = 0
            for i in range(len(U)):
                if (c >> i) & 1:
                    x ^= U[i]
            if all(dot(x, u) == 0 for u in U):
                cnt += 1
        assert cnt == 1 << len(rad)
    # hand cases: orthonormal family = label; hyperbolic plane = alt; isotropic vector = deg
    assert kind(rref([1, 2, 4]), 6) == "label"
    assert kind(rref([0b0011, 0b0110]), 4) == "alt"
    assert kind(rref([0b0011]), 4) == "deg"
    assert step_kind(rref([1]), rref([1, 2, 4]), 6) == "legal"
    assert step_kind(rref([1]), rref([1, 0b0110, 0b1100]), 6) == "alt"
    n = 0
    while n < 200:                  # orthonormal bases of random nondegenerate non-alternating subspaces
        U = rref([rnd.randrange(1, 1 << h) for _ in range(rnd.randrange(1, 9))])
        if kind(U, h) == "label":
            orthonormal(U)
            A = rref(orthonormal(U)[:rnd.randrange(0, len(U))])
            if kind(A, h) in ("zero", "label") and step_kind(A, U, h) == "legal":
                k, d = step_dirs(A, U, h)
                assert all(dot(x, a) == 0 for x in d for a in A) and rref(tuple(A) + tuple(d)) == U
            n += 1
    return "f2 selftest ok"


if __name__ == "__main__":
    print(selftest())
