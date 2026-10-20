"""x-transfer, 13 Oct 2026 (python3 -I, standard library).

P. E8 root lines (= the 120 odd-Y 4-qubit Paulis, x-dict section F): largest set of pairwise NON-orthogonal
   lines (= pairwise anticommuting Paulis), which bounds the copy-floor minrank from below over EVERY field;
   the Gram-rank bound 2n+1 (n = 4 here); and whether the 120 lines split into 15 mutually orthogonal octads
   (a cover by "commuting groups", an upper bound for the clique-cover decoder).
Q. The engine's "every f" rule versus single-copy quantum identities, at m = 1:
   - f = 1: K = SX is a rational combination of free matrices (so at one copy the kernel would cost nothing);
   - least f with K^(x f) outside the Q-span of {b^(x f) : b free};
   - least F such that Q[G_1] -> (+)_{f<=F} End((C^2)^(x f)) is injective (then 'correct for every f' is
     'correct for f = 0..F'; Burnside-Brauer).
"""
import itertools, sys
from fractions import Fraction

# ---------------- P ----------------
def e8_lines():
    roots = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 8; v[i] = si; v[j] = sj; roots.append(tuple(v))       # doubled coordinates
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            roots.append(tuple(signs))
    assert len(roots) == 240
    lines = {}
    for r in roots:
        k = max(r, tuple(-x for x in r))
        lines[k] = r
    return list(lines.values())

def ip(a, b): return sum(x * y for x, y in zip(a, b))   # 4 * true inner product

def max_clique(adj, n):
    best = []
    def expand(R, P):
        nonlocal best
        if not P:
            if len(R) > len(best):
                best = R[:]
            return
        if len(R) + len(P) <= len(best):
            return
        for v in sorted(P, key=lambda v: -len(adj[v] & P)):
            if len(R) + len(P) <= len(best):
                return
            expand(R + [v], P & adj[v])
            P = P - {v}
    expand([], set(range(n)))
    return best

def check_P():
    print("== P. E8 root lines: pairwise non-orthogonal sets, octad partitions")
    L = e8_lines()
    n = len(L)
    nonorth = [set(j for j in range(n) if j != i and ip(L[i], L[j]) != 0) for i in range(n)]
    orth = [set(j for j in range(n) if j != i and ip(L[i], L[j]) == 0) for i in range(n)]
    print(f"  {n} lines; each non-orthogonal to {sorted(set(len(s) for s in nonorth))}, orthogonal to {sorted(set(len(s) for s in orth))} others")
    C = max_clique(nonorth, n)
    print(f"  largest set of pairwise non-orthogonal lines: {len(C)} (Gram-rank bound for pairwise anticommuting 4-qubit Paulis: 2*4+1 = 9)")
    # verify the clique
    assert all(ip(L[a], L[b]) != 0 for a, b in itertools.combinations(C, 2))
    # octads: maximal sets of 8 mutually orthogonal lines
    octads = set()
    for i in range(n):
        def grow(S, cand):
            if len(S) == 8:
                octads.add(frozenset(S)); return
            for v in sorted(cand):
                if v > S[-1]:
                    grow(S + [v], cand & orth[v])
        grow([i], set(orth[i]))
    print(f"  octads (8 mutually orthogonal lines): {len(octads)}")
    # exact cover of the 120 lines by 15 octads
    by_line = {i: [o for o in octads if i in o] for i in range(n)}
    sol = []
    def cover(left, chosen):
        if not left:
            sol.append(list(chosen)); return True
        i = min(left, key=lambda v: sum(1 for o in by_line[v] if o <= left))
        for o in by_line[i]:
            if o <= left:
                chosen.append(o)
                if cover(left - o, chosen):
                    return True
                chosen.pop()
        return False
    found = cover(frozenset(range(n)), [])
    print(f"  partition of the 120 lines into 15 mutually orthogonal octads exists: {found}")
    if found:
        parts = sol[0]
        assert len(parts) == 15 and frozenset().union(*parts) == frozenset(range(n))
    return len(C)

# ---------------- Q ----------------
I2 = ((1, 0), (0, 1))
def kron(A, B):
    na, nb = len(A), len(B)
    return tuple(tuple(A[i // nb][j // nb] * B[i % nb][j % nb] for j in range(na * nb)) for i in range(na * nb))
def mul(A, B):
    n = len(A)
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)) for i in range(n))
def tpow(A, f):
    out = ((1,),)
    for _ in range(f):
        out = kron(A, out)
    return out
H = Fraction(1, 2)
SXm = ((complex(H, H), complex(H, -H)), (complex(H, -H), complex(H, H)))
Xm = ((0, 1), (1, 0)); Sm = ((1, 0), (0, 1j)); iI = ((1j, 0), (0, 1j))

def closure(gens):
    e = ((1 + 0j, 0j), (0j, 1 + 0j))
    seen = {e}; stack = [e]
    while stack:
        M = stack.pop()
        for g in gens:
            N = mul(g, M)
            if N not in seen:
                seen.add(N); stack.append(N)
    return list(seen)

def as_int_vec(M, f):
    s = 1 << f       # entries of an f-fold power of a 1-qubit element of G_1 have denominators dividing 2^f
    out = []
    for row in M:
        for v in row:
            v = complex(v)
            re, im = v.real * s, v.imag * s
            assert re == int(re) and im == int(im)
            out += [int(re), int(im)]
    return out

def rank_mod(rows, p):
    rows = [[x % p for x in r] for r in rows]
    rk = 0; ncol = len(rows[0]) if rows else 0; col = 0
    piv_rows = []
    for c in range(ncol):
        pr = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
        if pr is None:
            continue
        rows[rk], rows[pr] = rows[pr], rows[rk]
        inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [x * inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                t = rows[i][c]
                rows[i] = [(a - t * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
        if rk == len(rows):
            break
    return rk

def rank_exact(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        pr = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if pr is None:
            continue
        rows[rk], rows[pr] = rows[pr], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                t = rows[i][c] / rows[rk][c]
                rows[i] = [a - t * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
        if rk == len(rows):
            break
    return rk

def check_Q(fmax=5):
    print("== Q. 'every f' versus single-copy identities, m = 1")
    B1 = closure([Xm, Sm, iI])
    G1 = closure([Xm, Sm, iI, SXm])
    print(f"  |B_1| = {len(B1)}, |G_1| = {len(G1)}")
    for f in range(1, fmax + 1):
        span_rows = [as_int_vec(tpow(b, f), f) for b in B1]
        k = as_int_vec(tpow(SXm, f), f)
        r0 = rank_exact(span_rows)
        r1 = rank_exact(span_rows + [k])
        method = "exact over Q"
        print(f"  f={f}: dim_Q span{{b^(x f)}} = {r0}; with K^(x f) added: {r1}; K^(x f) in the free span: {r1 == r0}  [{method}]")
    # injectivity of Q[G_1] -> (+)_{f<=F} End
    for F in range(1, 9):
        rows = []
        for g in G1:
            v = []
            for f in range(F + 1):
                w = as_int_vec(tpow(g, f), f)
                # bring all levels to the common denominator 2^F
                v += [x << (F - f) for x in w]
            rows.append(v)
        rk = rank_mod(rows, 1000003)
        rkq = rank_exact(rows) if F in (4, 5) else None
        print(f"  F={F}: rank of the {len(G1)} vectors (g^(x f))_(f=0..F) mod p = {rk}" + (f", exact over Q = {rkq}" if rkq is not None else "") + ("  -> injective (rank mod p = 96 implies rank over Q = 96): 'every f' = 'f <= F'" if rk == len(G1) else ""))
        if rk == len(G1):
            break
    # distinct trace values (Burnside-Brauer predicts injectivity by F = (number of distinct values) - 1)
    tr = set()
    for g in G1:
        t = complex(g[0][0]) + complex(g[1][1])
        tr.add((round(t.real, 9), round(t.imag, 9)))
    print(f"  distinct trace values of G_1 in its defining representation: {len(tr)} (Burnside-Brauer bound F <= {len(tr) - 1})")

if __name__ == "__main__":
    check_P()
    check_Q()
    print("DONE")
