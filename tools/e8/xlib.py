"""xlib.py -- sums as bit masks, the order-64 group acting on them, orbits, pools.
A sum (symbol) is (mask, neg): mask = its sources (bit T), neg = the sources with coefficient -1; normalised so
that the lowest source has +1 (a sum is taken up to a global sign).  Standard library only."""
import sys
sys.dont_write_bytecode = True
from lib import *          # lib.py: e8, P, V, IDX, G, Fr, need, par
import gram, groups, hio

def low(m): return m & -m
def bits(m):
    out = []
    while m:
        b = m & -m; out.append(b.bit_length() - 1); m ^= b
    return out
def norm(m, n):
    if n & low(m): n ^= m
    return (m, n)
def sym_of(content):
    m = n = 0
    for T, x in content.items():
        m |= 1 << T
        if x < 0: n |= 1 << T
        assert abs(x) == 1
    return norm(m, n)
def content_of(s):
    m, n = s
    return {T: Fr(-1 if n >> T & 1 else 1) for T in bits(m)}

def group_elements(name="F8uu"):
    gens = [gram.lift(f) for f in groups.group(name)]
    def nz(p, e):
        if e[0] < 0: e = tuple(-x for x in e)
        return (tuple(p), tuple(e))
    gens = [nz(p, e) for p, e in gens]
    ident = (tuple(range(V)), tuple([1] * V))
    seen = {ident}; todo = [ident]
    while todo:
        p1, e1 = todo.pop()
        for p2, e2 in gens:
            g = nz([p2[p1[T]] for T in range(V)], [e1[T] * e2[p1[T]] for T in range(V)])
            if g not in seen: seen.add(g); todo.append(g)
    return sorted(seen)

class Act:
    """the group acting on sums and targets"""
    def __init__(self, name="F8uu"):
        self.els = group_elements(name)
        self.n = len(self.els)
        self.cache = {}
    def img(self, gi, s):
        p, e = self.els[gi]; m, n = s; m2 = n2 = 0; mm = m
        while mm:
            b = mm & -mm; T = b.bit_length() - 1; mm ^= b
            t2 = p[T]; m2 |= 1 << t2
            if (e[T] < 0) != bool(n >> T & 1): n2 |= 1 << t2
        return norm(m2, n2)
    def images(self, s):
        """list of the n images of s (index = group element)"""
        r = self.cache.get(s)
        if r is None:
            r = [self.img(gi, s) for gi in range(self.n)]
            self.cache[s] = r
        return r
    def orbit(self, s):
        return set(self.images(s))

NM = []; NN = []
for S in range(V):
    m = n = 0
    for T in range(V):
        if T != S and G[S][T]:
            m |= 1 << T
            if G[S][T] < 0: n |= 1 << T
    NM.append(m); NN.append(n)
def match(S, s):
    """+1 / -1 if target S can take sum s with coefficient +1/2 / -1/2, else 0"""
    m, n = s
    if m & ~NM[S]: return 0
    d = n ^ (NN[S] & m)
    if d == 0: return 1
    if d == m: return -1
    return 0
def targets_of(s):
    return [S for S in range(V) if match(S, s)]

class H: pass
def load_h(path):
    gd = hio.load_h(path)
    h = H(); h.content, h.frame, h.kids, h.need, h.twin = gd["content"], gd["frame"], gd["kids"], gd["need"], dict(gd.get("twin", {}))
    h.need = [dict(n) for n in h.need]
    return h
def needed_set(h):
    needed = set(); todo = [u for S in range(V) for u in h.need[S]]
    while todo:
        u = todo.pop()
        if u in needed: continue
        needed.add(u)
        if h.kids[u] is not None: todo.extend([h.kids[u][0], h.kids[u][1]])
    return needed
def proxy(h):
    """deliveries, ordinary sums, in-place pairs, sums with a single source as a part (ordinary)"""
    nd = needed_set(h)
    deliv = sum(len(x) for x in h.need)
    sums = [u for u in nd if h.kids[u] is not None]
    pr = set()
    for u in sums:
        t = h.twin.get(u)
        if t is not None and t in nd and h.kids[t][:2] == h.kids[u][:2]: pr.add(u)
    ordi = [u for u in sums if u not in pr]
    osrc = sum(1 for u in ordi if h.kids[u][0] < V or h.kids[u][1] < V)
    return dict(deliveries=deliv, sums=len(sums), pairs=len(pr) // 2, ordinary=len(ordi), ordinary_with_source=osrc, P=deliv + len(ordi))
