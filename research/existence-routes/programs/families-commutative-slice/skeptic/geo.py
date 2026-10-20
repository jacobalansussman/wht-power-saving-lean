# Own geometry for the skeptic check (written from scratch, nothing imported from ex2-families).
# Points of F_2^m are ints 1..2^m-1.  A "form support" is a set S of nonzero points whose indicator
# is a function of algebraic degree <= 2 vanishing at 0; its trace is #(basis vectors in S) mod 2.
from itertools import combinations


def anf_degree(ind, m):
    """algebraic degree of a boolean function given as a list of 2^m values"""
    a = list(ind)
    n = 1 << m
    for i in range(m):
        b = 1 << i
        for x in range(n):
            if x & b:
                a[x] ^= a[x ^ b]
    deg = 0
    for x in range(n):
        if a[x]:
            deg = max(deg, bin(x).count("1"))
    return deg


def span(vecs):
    s = {0}
    for v in vecs:
        s |= {x ^ v for x in s}
    return frozenset(s)


def subspaces(m, d):
    out = set()
    for basis in combinations(range(1, 1 << m), d):
        S = span(basis)
        if len(S) == (1 << d):
            out.add(S)
    return out


def flats_avoiding_0(m, d):
    out = set()
    for V in subspaces(m, d):
        for a in range(1, 1 << m):
            if a in V:
                continue
            out.add(frozenset(a ^ v for v in V))
    return sorted(out, key=lambda s: sorted(s))


def is_form_support(S, m):
    ind = [1 if x in S else 0 for x in range(1 << m)]
    return ind[0] == 0 and anf_degree(ind, m) <= 2


def trace(S, m):
    return sum(1 for j in range(m) if (1 << j) in S) % 2


def even_odd_flats(m, d):
    ev, od = [], []
    for S in flats_avoiding_0(m, d):
        assert is_form_support(S, m), "flat of dim %d at block %d not quadratic" % (d, m)
        (od if trace(S, m) else ev).append(S)
    return ev, od


def all_form_supports(m):
    """every binary quadratic form with Q(0)=0 as (support, trace); only for m <= 5"""
    pairs = [(j, k) for j in range(m) for k in range(j, m)]
    out = []
    for mask in range(1 << len(pairs)):
        q = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
        S = []
        for z in range(1, 1 << m):
            v = 0
            for (j, k) in q:
                v ^= (z >> j & 1) & (z >> k & 1)
            if v:
                S.append(z)
        tr = sum(1 for (j, k) in q if j == k) % 2
        out.append((frozenset(S), tr))
    return out
