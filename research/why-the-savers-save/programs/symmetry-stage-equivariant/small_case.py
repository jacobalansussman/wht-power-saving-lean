#!/usr/bin/env python3
"""small_case.py (d-I3-stage-equivariant): the hand-over question settled exactly on small analogues of the hub's
compressed-window cover (standard library only; python3 -I small_case.py).

Analogue: local F_2^h (h = 3 or 4), exterior blocks of size h-1, global m = h + (h-1)*S, dot form; ports = odd-weight
labels q with a zero coordinate (the hub's admissibility rule); routes rho_{s,q} built exactly as the hub's
INHERITED-COMPLEX-GEOMETRY (chart U_q, exchange q^perp with block s); cover = all of O(m,2), invocation d at stage s
processes port q's piece of class d*g_{s,q}.

 (A) h = 3, m = 5 (stages 0, 1): EVERY pair of routes (g_0, g_1) in O(5,2)^2 (720^2 pairs) is tested against the two
     boundary identities  g_0^-1(q^perp + E_0) = g_1^-1(E_1)  and  g_0^-1(U + E_0) = g_1^-1(<q> + E_1)  for each port.
     For every admissible pair, A = g_0 g_1^-1 maps E_1 onto q^perp; the sets of admissible A for different ports are
     disjoint.  So no route choice at all lets two ports of one invocation continue together.
 (B) h = 4, m = 7: the same with g_0 = I and g_1 over all of O(7,2) (1,451,520 elements).
 (C) h = 3, m = 7 (stages 0, 1, 2), the hub's own routes, full cover O(7,2) as invocations: symbolic replay of a
     stage-1 total over two ports pre-loaded into a third port's target register; at stage 2 the host's invocation
     never holds the two pieces the total carries (all invocations, all port triples).
"""
import itertools, sys, time

def dot(x, y): return bin(x & y).count('1') & 1
def rank(vs):
    b = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in b: v ^= b[h]
            else: b[h] = v; break
    return len(b)
def span_eq(A, B): return rank(A) == rank(B) == rank(list(A) + list(B))
def apply(cols, x):
    y = 0; i = 0
    while x:
        if x & 1: y ^= cols[i]
        x >>= 1; i += 1
    return y
def compose(A, B): return [apply(A, c) for c in B]
def inverse_orth(cols):          # orthogonal: inverse = transpose
    n = len(cols)
    return [sum(((cols[i] >> j) & 1) << i for i in range(n)) for j in range(n)]

def iter_orthogonal(m):
    """generate all ordered orthonormal bases of F_2^m (dot form) = O(m,2), as column tuples (streamed)"""
    unit = [x for x in range(1, 1 << m) if bin(x).count('1') & 1]
    def rec(cols):
        if len(cols) == m:
            yield tuple(cols); return
        for x in unit:
            if all(dot(x, c) == 0 for c in cols) and rank(cols + [x]) == len(cols) + 1:
                yield from rec(cols + [x])
    yield from rec([])
def orthogonal_group(m): return list(iter_orthogonal(m))

class Geometry:
    def __init__(self, h, S):
        self.h, self.bl, self.S = h, h - 1, S
        self.m = h + (h - 1) * S
        self.U = [1 << i for i in range(h)]
        self.ports = [q for q in range(1, 1 << h) if bin(q).count('1') & 1 and q != (1 << h) - 1]
    def block(self, s): return [1 << (self.h + self.bl * (s - 1) + k) for k in range(self.bl)]
    def E(self, s): return [v for j in range(1, s + 1) for v in self.block(j)]
    def perp(self, q):
        return [x for x in range(1, 1 << self.h) if dot(x, q) == 0]
    def chart(self, q):
        h = self.h
        z = next(i for i in range(h) if not (q >> i) & 1)
        a0 = q ^ (1 << z)
        order = list(range(h)); order[0], order[z] = order[z], order[0]
        B = [(1 << i) ^ (a0 if (a0 >> i) & 1 else 0) for i in order]
        assert B[0] == q and all(dot(B[i], B[j]) == int(i == j) for i in range(h) for j in range(h))
        return B
    def rho(self, s, q):
        m, h = self.m, self.h
        if s == 0: return [1 << i for i in range(m)]
        B = self.chart(q); bs = self.block(s)
        img = {B[0]: B[0]}
        for k in range(1, h): img[B[k]] = bs[k - 1]; img[bs[k - 1]] = B[k]
        cols = []
        for i in range(m):
            e = 1 << i
            if i < h:
                y = 0
                for k in range(h):
                    if dot(e, B[k]): y ^= img[B[k]]
                cols.append(y)
            else: cols.append(img.get(e, e))
        return cols
    def admissible(self, g0, g1, s, q):
        """boundary identities for port q between stage s (route g0) and s+1 (route g1), frames in class coordinates"""
        i0, i1 = inverse_orth(list(g0)), inverse_orth(list(g1))
        K = [apply(i0, v) for v in self.perp(q) + self.E(s)]
        Z = [apply(i1, v) for v in self.E(s + 1)]
        Hf = [apply(i0, v) for v in self.U + self.E(s)]
        P = [apply(i1, v) for v in [q] + self.E(s + 1)]
        return span_eq(K, Z) and span_eq(Hf, P)

def part_A():
    G = Geometry(3, 1); m = G.m
    t0 = time.time(); O = orthogonal_group(m)
    print('(A) h=3, m=%d, ports %s; |O(%d,2)| = %d (%.1fs)' % (m, G.ports, m, len(O), time.time() - t0))
    assert len(O) == 720
    Aset = {q: set() for q in G.ports}; npairs = {q: 0 for q in G.ports}
    hub_ok = {q: G.admissible(G.rho(0, q), G.rho(1, q), 0, q) for q in G.ports}
    for g0 in O:
        for g1 in O:
            for q in G.ports:
                if G.admissible(g0, g1, 0, q):
                    A = tuple(compose(list(g0), inverse_orth(list(g1))))
                    assert span_eq([apply(A, v) for v in G.E(1)], G.perp(q) + G.E(0))
                    Aset[q].add(A); npairs[q] += 1
    print('    hub routes admissible: %s' % hub_ok)
    print('    admissible (g_0, g_1) pairs per port: %s; distinct continuation factors A per port: %s'
          % (npairs, {q: len(v) for q, v in Aset.items()}))
    inter = max(len(Aset[a] & Aset[b]) for a, b in itertools.combinations(G.ports, 2))
    print('    every admissible A maps E_1 onto q^perp: True; largest overlap of the A-sets of two ports: %d' % inter)
    assert inter == 0 and all(hub_ok.values())

def part_B():
    """g_0 = I; for each g_1 in O(7,2): the boundary identities say g_1^-1(E_1) = q^perp (inside U) and
    g_1^-1(<q> + E_1) = U.  So g_1 determines q (perp of g_1^-1(E_1) inside U): each continuation factor serves at most
    one label.  Streamed over the whole group; counts the admissible g_1 per port and checks the direct test."""
    G = Geometry(4, 1); m = G.m; h = G.h
    t0 = time.time(); n = 0; per = {q: 0 for q in G.ports}; multi = 0
    Umask = (1 << h) - 1
    for g1 in iter_orthogonal(m):
        n += 1
        i1 = inverse_orth(list(g1))
        V = [apply(i1, v) for v in G.E(1)]
        if any(v & ~Umask for v in V): continue
        hits = [q for q in G.ports if G.admissible(tuple(1 << i for i in range(m)), g1, 0, q)]
        if len(hits) > 1: multi += 1
        for q in hits:
            assert span_eq(V, G.perp(q))
            per[q] += 1
    print('(B) h=4, m=%d, %d ports; streamed |O(%d,2)| = %d; admissible g_1 per port (g_0 = I): %s; '
          'g_1 admissible for two or more ports: %d (%.1fs)' % (m, len(G.ports), m, n, sorted(set(per.values())), multi, time.time() - t0))
    assert n == 1451520 and multi == 0 and all(per.values())

def random_orthogonal(m, rng):
    while True:
        cols = []
        unit = [x for x in range(1, 1 << m) if bin(x).count('1') & 1]
        while len(cols) < m:
            c = [x for x in unit if all(dot(x, y) == 0 for y in cols) and rank(cols + [x]) == len(cols) + 1]
            if not c: break
            cols.append(rng.choice(c))
        if len(cols) == m: return tuple(cols)

def part_C(sample=None):
    """sample=None: the whole cover O(7,2) (slow in pure Python); sample=N: d = I plus N random invocations.
    (Left multiplication by d is a bijection of the cover that commutes with every route, so d = I already decides it.)"""
    import random
    G = Geometry(3, 2); m = G.m
    t0 = time.time()
    R = {(s, q): G.rho(s, q) for s in range(3) for q in G.ports}
    bad = 0; checked = 0; s = 1; ncover = 0
    rng = random.Random(1013)
    cover = iter_orthogonal(m) if sample is None else [tuple(1 << i for i in range(m))] + [random_orthogonal(m, rng) for _ in range(sample)]
    for d in cover:
        ncover += 1
        for host in G.ports:
            others = [q for q in G.ports if q != host]
            # pieces carried by a stage-1 total over the other ports, built in invocation d: class d*rho_{1,q}
            carried = {(q, tuple(compose(list(d), R[s, q]))) for q in others}
            # the host's stage-2 invocation d' = d * rho_{1,host} * rho_{2,host}; the pieces it processes: d'*rho_{2,q}
            d2 = compose(compose(list(d), R[s, host]), R[s + 1, host])
            present = {(q, tuple(compose(d2, R[s + 1, q]))) for q in G.ports}
            assert (host, tuple(compose(list(d), R[s, host]))) in present        # the host's own piece continues
            if carried & present: bad += 1
            checked += 1
    print('(C) h=3, m=%d, stages 0..2, hub routes, %s: %d invocations, (invocation, host) cases %d; '
          'cases where a carried piece is present in the host continuation: %d (%.1fs)'
          % (m, 'whole cover O(%d,2)' % m if sample is None else 'd = I plus %d random invocations in O(%d,2)' % (sample, m),
             ncover, checked, bad, time.time() - t0))
    assert bad == 0 and (sample is not None or ncover == 1451520)

if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'C', 'B']
    for w in which: {'A': part_A, 'B': part_B, 'C': part_C, 'Cs': lambda: part_C(2000)}[w]()
    print('PASS')
