"""sk_census.py -- skeptic v-tb4-nonlinear, independent recount of widened-game shortest routes to the kernel.

Independent model (no complex arithmetic, no AGL orbits). K = Phi * H * Phi' with H the Sylvester +-1 matrix and
Phi, Phi' diagonal Z4 (rows of K: (-i)^{wt(x^y)} = (-i)^{wt x} (-i)^{wt y} (-1)^{x.y}). Left monomials do not change
a coset B g; the common right diagonal Phi' commutes with every backward step (left multiplications), so the coset
chains of K and of H correspond one to one. A route row is c * chi_p restricted to a support S; a row is stored as
(S, N): S = support bitmask over the 2^m points, N = points of S where the row has the opposite sign to its value at
min(S). Backward step (D^-1 b^-1 g): perfect matching of rows into pairs with EQUAL support whose ratio
(N1 xor N2 inside S) is balanced; children = row 1 restricted to the two ratio classes. Frame = sorted multiset of
rows (row up to scalar; scalars of route rows always differ by Z4 since every step multiplies by e^{+-i pi/4} sqrt2).
Exhaustive over all perfect matchings of K's rows; memo keyed by exact bytes of the frame."""
import sys, time, itertools
from collections import Counter

m = int(sys.argv[1])
n = 1 << m
FULL = (1 << n) - 1
popc = lambda v: bin(v).count("1")


def norm(S, N):
    low = S & -S
    if N & low:
        N ^= S
    return (S, N)


def krows():
    rows = []
    for x in range(n):
        N = 0
        for y in range(n):
            if popc(x & y) & 1:
                N |= 1 << y
        rows.append(norm(FULL, N))
    return rows


def key(rows):
    return b"".join(S.to_bytes(4, "little") + N.to_bytes(4, "little") for (S, N) in sorted(rows))


def step_pair(r1, r2):
    S, N1 = r1
    _, N2 = r2
    R = (N1 ^ N2) & S
    P = S & ~R
    return norm(P, N1 & P), norm(R, N1 & R)


def compatible(r1, r2):
    S = r1[0]
    if S != r2[0]:
        return False
    R = (r1[1] ^ r2[1]) & S
    return 2 * popc(R) == popc(S)


def group_matchings(rows):
    """yield all perfect matchings (lists of index pairs) of rows into compatible pairs"""
    k = len(rows)
    used = [False] * k
    cur = []

    def rec():
        try:
            p = used.index(False)
        except ValueError:
            yield list(cur)
            return
        used[p] = True
        for q in range(p + 1, k):
            if not used[q] and compatible(rows[p], rows[q]):
                used[q] = True
                cur.append((p, q))
                yield from rec()
                cur.pop()
                used[q] = False
        used[p] = False
    yield from rec()


def children(rows):
    """distinct child frames (key -> rows)"""
    groups = {}
    for r in rows:
        groups.setdefault(r[0], []).append(r)
    glist = list(groups.values())
    if any(len(g) % 2 for g in glist):
        return {}
    per = []
    for g in glist:
        opts = []
        for M in group_matchings(g):
            ch = []
            for (p, q) in M:
                ch.extend(step_pair(g[p], g[q]))
            opts.append(ch)
        if not opts:
            return {}
        per.append(opts)
    out = {}
    for combo in itertools.product(*per):
        rows2 = [r for part in combo for r in part]
        out.setdefault(key(rows2), rows2)
    return out


memo = {}       # key -> number of routes from this frame down to depth 0
levels = {}     # depth -> set of keys with >0 routes


def nroutes(rows, d):
    if d == 0:
        assert sorted(S for S, _ in rows) == [1 << i for i in range(n)], "not monomial at depth 0"
        return 1
    k = key(rows)
    if k in memo:
        return memo[k]
    tot = 0
    for kk, r2 in children(rows).items():
        tot += nroutes(r2, d - 1)
    if d < m - 1:
        memo[k] = tot
    if tot:
        levels.setdefault(d, set()).add(k)
    return tot


t0 = time.time()
K = krows()
assert len(set(K)) == n
# top level: every perfect matching of K's rows, with no compatibility assumption (we check it)
allpairs_ok = all(compatible(K[p], K[q]) for p in range(n) for q in range(p + 1, n))
print("m", m, "all pairs of K rows compatible:", allpairs_ok, flush=True)
nmatch = 0
ncont = 0
ntrans_cont = 0
total = 0
trans_routes = 0
child_count = Counter()     # child key -> number of matchings giving it (continuing only)
child_routes = {}
for M in group_matchings(K):
    nmatch += 1
    ch = []
    for (p, q) in M:
        ch.extend(step_pair(K[p], K[q]))
    kk = key(ch)
    r = nroutes(ch, m - 1)
    if r:
        ncont += 1
        child_count[kk] += 1
        child_routes[kk] = r
        diffs = {p ^ q for (p, q) in M}      # rows of H indexed by x; K row order = x
        if len(diffs) == 1:
            ntrans_cont += 1
            trans_routes += r
    if nmatch % 200000 == 0:
        print(" ..", nmatch, ncont, round(time.time() - t0, 1), "memo", len(memo), flush=True)
print("perfect matchings of K rows:", nmatch, flush=True)
print("continuing first steps (matchings):", ncont, "distinct depth-(m-1) frames:", len(child_count),
      "max matchings per frame:", max(child_count.values()), flush=True)
print("continuing translation matchings:", ntrans_cont, "routes through them:", trans_routes, flush=True)
total = sum(child_routes.values())
print("TOTAL ROUTES (chains of distinct frames):", total, flush=True)
for d in sorted(levels):
    print("  |S_%d| = %d" % (d, len(levels[d])), flush=True)
print("t=%.1f memo=%d" % (time.time() - t0, len(memo)), flush=True)
memo.clear(); levels.clear(); child_count.clear(); child_routes.clear()


# ---- control: routes whose every frame has all row supports cosets of ONE subspace (the affine game's flags)
def dirspace(S):
    pts = [y for y in range(n) if S >> y & 1]
    return frozenset(a ^ b for a in pts for b in pts)


def uniform(rows):
    ds = {dirspace(S) for S, _ in rows}
    if len(ds) != 1:
        return False
    V = next(iter(ds))
    # V must be closed under xor (a subspace) and every support a coset of it
    return all((a ^ b) in V for a in V for b in V) and all(len(V) == bin(S).count("1") for S, _ in rows)


amemo = {}


def nroutes_aff(rows, d):
    if d == 0:
        return 1
    k = key(rows)
    if k in amemo:
        return amemo[k]
    tot = 0
    groups = {}
    for r in rows:
        groups.setdefault(r[0], []).append(r)
    glist = list(groups.values())
    if any(len(g) % 2 for g in glist):
        amemo[k] = 0
        return 0
    def gen_opts(g):
        for M in group_matchings(g):
            ch = []
            for (p, q) in M:
                ch.extend(step_pair(g[p], g[q]))
            yield ch
    if len(glist) == 1:          # stream (top level: 2M matchings, do not store them)
        combos = ((ch,) for ch in gen_opts(glist[0]))
    else:
        combos = itertools.product(*[list(gen_opts(g)) for g in glist])
    seen = set()
    for combo in combos:
        r2 = [r for part in combo for r in part]
        if not uniform(r2):
            continue
        kk = key(r2)
        if kk in seen:
            continue
        seen.add(kk)
        tot += nroutes_aff(r2, d - 1)
    amemo[k] = tot
    return tot


print("CONTROL affine-flag routes:", nroutes_aff(K, m), "(expected: number of complete flags of F_2^%d)" % m)
