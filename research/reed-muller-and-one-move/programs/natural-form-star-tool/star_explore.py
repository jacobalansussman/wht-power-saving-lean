"""star_explore.py <k list>: for one canonical target set, the star-lemma recursion (e3-slice6 7.2) on every even form
F: u(S) for S = D_F n (affine flat), dims 4 -> 0, starting from u(S) = min(|S|, s_G over forms G != F with S in D_G),
then u(D_F n A) <= floor((sum over the g flats B > A of dim+1 of u(D_F n B) - 2|D_F n A|) / (g - 2)).
A negative value is a contradiction (no program). Prints statistics. (Exploration; any hit is re-checked separately.)"""
import sys, struct, itertools, collections
import numpy as np
ks = [int(x) for x in sys.argv[1].split(",")]
raw = open("tables5.bin", "rb").read()
nf = struct.unpack_from("<i", raw, 0)[0]; o = 4
mask = np.frombuffer(raw, "<u4", nf, o).astype(np.int64); o += 4 * nf
wt = np.frombuffer(raw, "<i4", nf, o).astype(np.int64); o += 4 * nf
vid = np.frombuffer(raw, "<i4", nf, o).astype(np.int64); o += 4 * nf
o += 4; pk = np.frombuffer(raw, "<u4", 32, o).astype(np.int64)
cnt = np.zeros(31, dtype=np.int64)
for k in ks: cnt[k % 31] += 1
w = np.array([sum(cnt[k] for k in range(31) if pk[v] >> k & 1) for v in range(32)])
s = wt - w[vid]
assert s[1:].min() >= 0
# affine flats of F_2^5 as bitmasks over directions 1..31 (bit z-1); point 0 dropped
def pts_mask(P): return sum(1 << (z - 1) for z in P if z)
subs = {0: [frozenset([0])]}
allv = range(32)
def span(gens):
    S = {0}
    for g in gens: S |= {x ^ g for x in S}
    return frozenset(S)
for d in range(1, 6):
    seen = set()
    for gens in itertools.combinations(range(1, 32), d):
        V = span(gens)
        if len(V) == 1 << d: seen.add(V)
    subs[d] = list(seen)
flats = {}
for d in range(0, 5):
    L = set()
    for V in subs[d]:
        for a in range(32): L.add(frozenset(a ^ v for v in V))
    flats[d] = list(L)
print("flats per dim", {d: len(flats[d]) for d in flats})
up = {}
for d in range(0, 4):
    idx = {f: i for i, f in enumerate(flats[d + 1])}
    up[d] = [[idx[B] for B in flats[d + 1] if A <= B] for A in flats[d]]
fm = {d: np.array([pts_mask(A) for A in flats[d]], dtype=np.int64) for d in flats}
even = [F for F in range(1, nf) if vid[F] == 0]
Gm = mask[1:]; Gs = s[1:]
res = collections.Counter(); worst = []
for F in even:
    DF = int(mask[F]); n = bin(DF).count("1")
    u = {}
    for d in range(4, -1, -1):
        S = fm[d] & DF
        sz = np.array([bin(int(x)).count("1") for x in S])
        base = np.minimum(sz, 10 ** 6)
        for i, x in enumerate(S):
            x = int(x)
            if x == 0: base[i] = 0; continue
            ok = (Gm & x) == x
            ok[F - 1] = False
            if ok.any(): base[i] = min(base[i], int(Gs[ok].min()))
        u[d] = base.copy()
        if d < 4:
            g = (1 << (5 - d)) - 1
            for i in range(len(S)):
                tot = sum(int(u[d + 1][j]) for j in up[d][i])
                val = (tot - 2 * int(sz[i]))
                if val < 0:
                    res["contradiction"] += 1; worst.append((F, d, i, tot, int(sz[i]))); break
                u[d][i] = min(u[d][i], val // (g - 2))
            else:
                continue
            break
    # top: the 62 half-sets vs the whole (d = 5): star at d=4 is not defined (g = 1); report min margins at d = 3
print("even forms", len(even), dict(res), worst[:5])
