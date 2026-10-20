# v-iso-combi check v1 (skeptic of iso-combi). Small, pure Python, seconds, < 50 MB.
# (a) explicit graph isomorphism: hub "E8" labels (know = even overlap)  <->  E8 root lines (non-orthogonal)
# (b) explicit graph isomorphism: block-3 frames (opposite)  <->  singular points of the E8-mod-2 quadric (non-orthogonal)
#     (iso-combi only compared SRG parameters for (b) and cited triality from memory)
# (c) the order-3 separation identity of v-ms-fano: 2 pieces suffice outside 2-groups (limits the RM(1,k) picture)
# (d) Hirzebruch arithmetic for AG(2,q) dualised (t_q = q^2 + q lines of q points, n = q^2 points)
from itertools import combinations, product
import sys, time

def find_iso(n1, adj1, n2, adj2, tlimit=60.0):
    """Backtracking with forward checking on bitsets. adj*: list of int bitmasks. Returns list f or None."""
    if n1 != n2: return None
    full = (1 << n2) - 1
    non2 = [full & ~adj2[w] & ~(1 << w) for w in range(n2)]
    t0 = time.time()
    f = [-1] * n1
    def rec(cand, unm):
        if time.time() - t0 > tlimit: raise TimeoutError
        if not unm: return True
        u = min(unm, key=lambda x: bin(cand[x]).count("1"))
        c = cand[u]
        while c:
            w = (c & -c).bit_length() - 1; c &= c - 1
            new = dict(cand); ok = True
            for u2 in unm:
                if u2 == u: continue
                m = adj2[w] if (adj1[u] >> u2) & 1 else non2[w]
                new[u2] = cand[u2] & m & ~(1 << w)
                if not new[u2]: ok = False; break
            if not ok: continue
            f[u] = w
            if rec(new, [x for x in unm if x != u]): return True
            f[u] = -1
        return False
    cand = {u: full for u in range(n1)}
    sys.setrecursionlimit(10000)
    return f if rec(cand, list(range(n1))) else None

def check_iso(f, n, adj1, adj2):
    return sorted(f) == list(range(n)) and all(((adj1[a] >> b) & 1) == ((adj2[f[a]] >> f[b]) & 1) for a in range(n) for b in range(n) if a != b)

def to_adj(V, rel):
    idx = {v: i for i, v in enumerate(V)}
    adj = [0] * len(V)
    for i, a in enumerate(V):
        for j, b in enumerate(V):
            if i != j and rel(a, b): adj[i] |= 1 << j
    return adj

wt = lambda x: bin(x).count("1")
# (a)
labels = [x for x in range(1 << 9) if wt(x) in (3, 7)]
adjL = to_adj(labels, lambda x, y: wt(x & y) % 2 == 0)
roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product((2, -2), repeat=2):
        r = [0] * 8; r[i] = si; r[j] = sj; roots.append(tuple(r))
for s in product((1, -1), repeat=8):
    if sum(1 for t in s if t < 0) % 2 == 0: roots.append(s)
lines = sorted({max(r, tuple(-t for t in r)) for r in roots})
ip = lambda a, b: sum(x * y for x, y in zip(a, b)) // 4
adjR = to_adj(lines, lambda a, b: ip(a, b) != 0)
t = time.time(); f = find_iso(120, adjL, 120, adjR)
print("(a) labels(know) ~ E8 root lines(non-orthogonal): iso found:", f is not None, "verified:", f is not None and check_iso(f, 120, adjL, adjR), f"({time.time()-t:.1f}s)")
# (b)
def span(vs):
    S = {0}
    for v in vs: S |= {s ^ v for s in S}
    return frozenset(S)
def symp(a, b):
    return (wt((a & 7) & (b >> 3)) + wt((a >> 3) & (b & 7))) & 1
L = set()
for a, b, c in combinations(range(1, 64), 3):
    if symp(a, b) or symp(a, c) or symp(b, c): continue
    S = span([a, b, c])
    if len(S) == 8: L.add(S)
L = sorted(L, key=sorted)
adjF = to_adj(L, lambda A, B: len(A & B) == 1)          # opposite frames: trivial intersection
E = [x for x in range(1 << 9) if wt(x) % 2 == 0]
sing = [x for x in E if x and (wt(x) // 2) % 2 == 0]
adjS = to_adj(sing, lambda x, y: wt(x & y) % 2 == 1)    # non-orthogonal singular points
t = time.time()
try:
    f = find_iso(len(L), adjF, len(sing), adjS)
    print("(b) frames(opposite) ~ singular points of O+(8,2)(non-orthogonal): iso found:", f is not None,
          "verified:", f is not None and check_iso(f, len(L), adjF, adjS), f"({time.time()-t:.1f}s)")
except TimeoutError:
    print("(b) search timed out")
# also: distance-1 frames (14 per frame) should become a structure on the quadric side; report the image
# of 'frames at distance 1' as a relation between singular points: collinear (orthogonal) or not
if f is not None:
    inv = {f[i]: i for i in range(len(L))}
    d1 = [(i, j) for i in range(len(L)) for j in range(len(L)) if i < j and len(L[i] & L[j]) == 4]
    kinds = {}
    for i, j in d1:
        x, y = sing[f[i]], sing[f[j]]
        kinds[wt(x & y) % 2] = kinds.get(wt(x & y) % 2, 0) + 1
    print("    image of distance-1 frame pairs under the found map: B(x,y) values", kinds, "(0 = orthogonal/collinear)")
# (c) order-3 separation: K = H^{(x)3} (integer, unnormalised), sigma = cyclic shift of 3 address bits.
N = 8
H = [[(-1) ** wt(r & c) for c in range(N)] for r in range(N)]
def rot(x): return ((x << 1) & 7) | (x >> 2)
Sg = [[1 if r == rot(c) else 0 for c in range(N)] for r in range(N)]
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
def add(A, B, s=1): return [[A[i][j] + s * B[i][j] for j in range(N)] for i in range(N)]
I = [[int(i == j) for j in range(N)] for i in range(N)]
S2 = mm(Sg, Sg)
twoZ2 = mm(H, add(add(I, Sg, -1), S2))           # 2*Z2 = K (1 - s + s^2)
lhs = mm(twoZ2, add(I, Sg))                      # 2*(Z1 + Z2) = 2 Z2 (1 + s)  must equal 2K
print("(c) order-3 letter: 2 Z2 (1 + sigma) == 2 K:", lhs == [[2 * x for x in row] for row in H], "; sigma^3 = 1:", mm(S2, Sg) == I)
# with only Z^v (2-group) relabellings, two pieces with distinct relabellings are never frustrated:
print("    2-group case: pieces x1 + s1 Z^u x2, x1 + s2 Z^v x2 (u != v): any sign pair is +-chi_a on {u,v} -> not separable")
# (d)
for q in (3, 4, 5, 7, 8, 9):
    n, tq = q * q, q * q + q
    lhs = (0.75 * tq) if q == 3 else 0.0
    rhs = n + (q - 4) * tq if q >= 5 else n
    print(f"(d) AG(2,{q}): n={n}, t_{q}={tq}: t2 + 3/4 t3 = {lhs}  vs  n + sum (r-4) t_r = {rhs}:", "holds" if lhs >= rhs else "violated")
