"""toy.py (u-ours, wht43/copywall) -- exact small model of the copy wall on E8 and on the hub's 112-line puncture.

Model (ours, wht35/e8/s-copies, 'reads at 0'): ports = 120 E8 labels in F_2^9 (9-point model: triples and
complements of pairs); S knows T iff even overlap.  A total g_k reads the sources T with g_k[T] != 0 and stands at a
frame F_k (subspace of F_2^9) containing every label it reads; its copy costs dim F_k.  A decoder is a family of
totals with scatter coefficients r_k such that B[S][T] = sum_k r_k[S] g_k[T] has B[S][S] = 1 and B[S][T] = 0 whenever
S does not know T (S != T).  cst = sum_k dim F_k.

P112 = E8 minus one octad (8 mutually orthogonal roots).  In the Hamming realisation (E8 from the extended Hamming
code) the 8 coordinate lines +-2e_i are such an octad and the 112 sign lines are the rest, so P112 is the hub's
complex unit (#523 -> #799 lineage, 112 ports, copy cost 64), up to the E8 symmetry.

Parts: (a) the 8 'coordinate totals' of the removed octad fit P112 exactly at cst 64 (hub control);
(b) octad-count floor for P112 over ALL subspaces (our 55-argument transferred); (c) uniform-dimension table;
(d) linear (natural-family) totals below dimension 8.
usage: python3 -I toy.py [octad index, default 0]"""
import sys, itertools, collections, time
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), "../e8-labels"))
import e8
import numpy as np
from fractions import Fraction as Fr

V = e8.V; PORTS = e8.PORTS
GRAM = np.array(e8.GRAM, dtype=np.int64)
ORTH = (GRAM == 0); KNOW = (np.abs(GRAM) == 1)
def par(x): return bin(x).count("1") & 1

def cliques(adj, size, universe):
    out = []
    def rec(cur, cand):
        if len(cur) == size:
            out.append(tuple(cur)); return
        c = cand
        while c:
            low = c & -c; i = low.bit_length() - 1; c ^= low
            if bin(cand >> i).count("1") + len(cur) < size: break
            rec(cur + [i], cand & adj[i] & ~((low << 1) - 1))
    rec([], universe)
    return out

ORTHMASK = [sum(1 << T for T in range(V) if ORTH[S, T]) for S in range(V)]
OCT = cliques(ORTHMASK, 8, (1 << V) - 1)
assert len(OCT) == 2025
oi = int(sys.argv[1]) if len(sys.argv) > 1 else 0
O = OCT[oi]
P = [T for T in range(V) if T not in O]
Pset = set(P)
print("octad removed:", O, [e8.NAMES[T] for T in O] if hasattr(e8, "NAMES") else "")
print("|P| =", len(P))

def span_dim(labels): return len(e8.echelon(list(labels)))

# (a) coordinate totals of the removed octad
ok = True
cst = 0
for o in O:
    reads = [T for T in P if GRAM[o, T] != 0]
    d = span_dim(PORTS[T] for T in reads)
    cst += d
for S in P:
    for T in P:
        b = sum(Fr(int(GRAM[o, S]) * int(GRAM[o, T]), 4) for o in O)
        want = Fr(int(GRAM[S, T]), 2)
        if b != want: ok = False
print("(a) 8 coordinate totals of the removed octad: B = Gram/2 on P exactly: %s; reads per total %s; cst = %d"
      % (ok, sorted(collections.Counter(int((GRAM[o, P] != 0).sum()) for o in O).items()), cst))
# same totals on the full E8 (each must also read its own octad port): control 72
cst120 = sum(span_dim(PORTS[T] for T in range(V) if GRAM[o, T] != 0) for o in O)
print("    the same totals on all 120 ports (they then read o itself): cst = %d (control: 72)" % cst120)

# octads inside P
OCTP = [oc for oc in OCT if Pset.issuperset(oc)]
print("octads inside P:", len(OCTP), " (each port of P in %s of them)" % sorted(collections.Counter(T for oc in OCTP for T in oc).values().__iter__().__class__ and collections.Counter(collections.Counter(T for oc in OCTP for T in oc).values()).items()))
inc = np.zeros((V, len(OCTP)), dtype=np.float32)
for j, oc in enumerate(OCTP):
    for T in oc: inc[T, j] = 1
Pmask = np.zeros(V, dtype=bool); Pmask[P] = True

# (b) all subspaces F of F_2^9 as C-perp; ports of P in F; octads (inside P) met
P0 = np.array([[par(PORTS[T] & c) == 0 and Pmask[T] for T in range(V)] for c in range(512)], dtype=bool)
P0pk = np.packbits(P0, axis=1)
POP = np.array([bin(i).count("1") for i in range(256)], dtype=np.int64)
def subspaces(k, n=9):
    for piv in itertools.combinations(range(n - 1, -1, -1), k):
        free = [[j for j in range(p) if j not in piv] for p in piv]
        f = sum(len(x) for x in free)
        idx = np.arange(1 << f, dtype=np.int64); cols = []; sh = 0
        for p, fr in zip(piv, free):
            col = np.full(1 << f, 1 << p, dtype=np.int64)
            for j in fr:
                col |= ((idx >> sh) & 1) << j; sh += 1
            cols.append(col)
        yield np.stack(cols, axis=1)

MAXHIT = {}; MAXPORT = {}; ROWS = {}
for d in (8, 7, 6, 5, 4, 3, 2, 1):
    t0 = time.time(); k = 9 - d; seen = set(); pmax = 0
    for bas in subspaces(k):
        for a in range(0, len(bas), 100000):
            b = bas[a:a + 100000]
            m = P0pk[b[:, 0]]
            for i in range(1, k): m = m & P0pk[b[:, i]]
            cnt = POP[m].sum(1); pmax = max(pmax, int(cnt.max()))
            for row in np.unique(m, axis=0): seen.add(row.tobytes())
    keys = list(seen)
    M = np.unpackbits(np.frombuffer(b"".join(keys), dtype=np.uint8).reshape(len(keys), 15), axis=1)[:, :V].astype(bool)
    hits = np.zeros(len(keys), dtype=np.int64)
    for a in range(0, len(keys), 2000):
        hits[a:a + 2000] = ((M[a:a + 2000].astype(np.float32) @ inc) > 0.5).sum(1)
    incb = inc > 0.5
    for j in list(np.argsort(-hits)[:200]):
        assert int(incb[M[j]].any(0).sum()) == hits[j]
    MAXHIT[d] = int(hits.max()); MAXPORT[d] = pmax; ROWS[d] = M
    print("(b) d = %d: %d distinct port sets, most ports %d, most octads met %d of %d (%.1f per dim) [%.0f s]"
          % (d, len(keys), pmax, MAXHIT[d], len(OCTP), MAXHIT[d] / d, time.time() - t0), flush=True)
MAXHIT[9] = len(OCTP); MAXPORT[9] = len(P)
NEED = 8 * len(OCTP)
dp = [0] + [-1] * 200
for c in range(1, 201):
    dp[c] = max(dp[c - d] + MAXHIT[d] for d in MAXHIT if d <= c)
cmin = min(c for c in range(201) if dp[c] >= NEED)
print("(b) octad-count floor on P: need %d octad-hits; LOWER BOUND cst >= %d (most hits at %d: %d)" % (NEED, cmin, cmin - 1, dp[cmin - 1]))

# (c) uniform dimension table: Hall (octads), reads (>= 2|P| - cst), rank (>= 8)
v = len(P)
print("(c) uniform designs on P (all totals at one dimension d):")
for d in range(1, 10):
    kh = -(-NEED // MAXHIT[d]) if MAXHIT[d] else 999
    kr = next(K for K in range(1, 999) if MAXPORT[d] * K >= 2 * v - d * K)
    kk = max(kh, kr, 8)
    print("    d = %d: ports/frame <= %3d, octads met <= %4d; K >= %d (Hall), %d (reads), 8 (rank): cost >= %d; K with cost < 64: %s"
          % (d, MAXPORT[d], MAXHIT[d], kh, kr, d * kk, [K for K in range(kk, 80) if d * K < 64]))

# (d) linear totals: see toy_d.py (the float rank loop that stood here was too slow).
