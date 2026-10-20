"""fibcount.py (v-t2-pair): number of H-orbits of nonzero binary quadratic forms (the fibres) in the lifted game
(Singer-15 on the first 4 bits, identity on bit 5) and at blocks 4, 5. Burnside over the cyclic group; pure Python."""
import itertools
def field(m, poly):
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> m: a ^= poly
        return r
    return mul
def count(m, poly, extra):
    mul = field(m, poly); q = 1 << m; n = q - 1; mm = m + extra
    pairs = [(i, j) for i in range(mm) for j in range(i, mm)]
    def bits(x): return [(x >> i) & 1 for i in range(mm)]
    pts = list(range(1 << mm))
    def act(x):  # sigma on the low m bits, identity on the rest
        lo = x & (q - 1); hi = x >> m
        return mul(lo, 2) | (hi << m)
    perm = [act(x) for x in pts]
    tabs = set(); orbits = 0; seen = set()
    for qv in itertools.product((0, 1), repeat=len(pairs)):
        if not any(qv): continue
        T = tuple(sum(c & bits(x)[i] & bits(x)[j] for c, (i, j) in zip(qv, pairs)) & 1 for x in pts)
        if T in seen: continue
        orbits += 1; cur = T
        for _ in range(n):
            seen.add(cur); cur = tuple(cur[perm[x]] for x in pts)
    return orbits
print("block 4 (Singer-15):", count(4, 0b10011, 0))
print("lifted block 5 (Singer-15 on 4 bits + fixed bit):", count(4, 0b10011, 1))
