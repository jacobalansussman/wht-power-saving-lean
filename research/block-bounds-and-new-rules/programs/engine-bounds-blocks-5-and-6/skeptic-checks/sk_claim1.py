# v-pf-induction-compute: explicit test of PROOF.md Theorem 3 Claim 1 (abar^(k)_r = sum_s N_r(s) abar^(k-1)_s)
# on an explicitly BUILT projected family, with my own symplectic identifications iota_p.
# Base: a random coordinate family at block 2 (random subsets S(U) of {0..19}, no symmetry) and a random
# subspace family over GF(10007) at block 2 (dims 3..5 in a 6-dim space; delta by ranks mod p; fine since the
# identity is linear in delta, the field does not matter). Project to block 3 explicitly, compute abar^(3)_r over
# all ordered pairs, compare exactly. Then repeat block 3 -> block 4 on a sample-free exact route for the
# coordinate family (all 2295^2 pairs is too many in Python, so block 3 -> 4 is skipped).
import itertools, random
from fractions import Fraction as Fr
def form(m, a, b):
    M = (1 << m) - 1
    return (bin(a & M & (b >> m)).count('1') + bin(b & M & (a >> m)).count('1')) & 1
def span(gens):
    S = {0}
    for g in gens:
        if g not in S: S |= {s ^ g for s in S}
    return frozenset(S)
def lagrangians(m):
    out = set()
    for B in itertools.combinations(range(1, 1 << (2 * m)), m):
        if all(form(m, x, y) == 0 for x in B for y in B):
            S = span(B)
            if len(S) == 2 ** m: out.add(S)
    return sorted(out, key=lambda S: sorted(S))
def dist(m, A, B): return m - (len(A & B).bit_length() - 1)
def iota(k, p):
    """a symplectic iso from p^perp/p (block k) to standard F_2^(2(k-1)); returns w -> int."""
    n = 2 * k
    q = next(v for v in range(1, 1 << n) if form(k, p, v))
    W = [v for v in range(1, 1 << n) if form(k, v, p) == 0 and form(k, v, q) == 0]   # p^perp cap q^perp
    E, F = [], []
    def red(v):
        for e, f in zip(E, F):
            if form(k, v, f): v ^= e
            if form(k, v, e): v ^= f
        return v
    pool = list(W)
    while len(E) < k - 1:
        e = next(red(v) for v in pool if red(v))
        f = next(red(v) for v in pool if red(v) and form(k, red(v), e))
        f = red(f); E.append(e); F.append(f)
    def phi(w):
        assert form(k, w, p) == 0
        w ^= p if form(k, w, q) else 0           # component in p^perp cap q^perp (kills p)
        x = sum(form(k, w, F[i]) << i for i in range(k - 1))
        z = sum(form(k, w, E[i]) << i for i in range(k - 1))
        return x | (z << (k - 1))
    return phi
def proj_index(k, phi, p, U, index_lower):
    A = span([v for v in U if form(k, v, p) == 0] + [p])
    img = frozenset(phi(v) for v in A)
    assert len(img) == 2 ** (k - 1)
    return index_lower[img]
def rank_mod(rows, P=10007):
    rows = [list(r) for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % P), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], P - 2, P)
        rows[rk] = [x * inv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % P:
                t = rows[i][c]; rows[i] = [(a - t * b) % P for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def closed(m, r):
    d = {}
    for s, c in ((r, (2**(m-r) - 1) + (4**m - 2**(m+r))),
                 (r-1, 2*(2**m - 2**(m-r)) + (2**r - 1) * 2**(r-1) * 2**(m-r)),
                 (r-2, (2**r - 1) * (2**(r-1) - 1) * 2**(m-r))):
        if c: d[s] = d.get(s, 0) + c
    return d
rng = random.Random(4242)
L2 = lagrangians(2); L3 = lagrangians(3)
print(f"Lagrangians: block 2 {len(L2)}, block 3 {len(L3)}")
idx2 = {S: i for i, S in enumerate(L2)}
phis = {p: iota(3, p) for p in range(1, 64)}
proj = {(u, p): proj_index(3, phis[p], p, U, idx2) for u, U in enumerate(L3) for p in range(1, 64)}
# base families: delta tables on block-2 frames
famS = [frozenset(x for x in range(20) if rng.random() < 0.45) for _ in L2]
dC = [[len(famS[a] - famS[b]) for b in range(len(L2))] for a in range(len(L2))]
P = 10007; D = 6
famV = [[[rng.randrange(P) for _ in range(D)] for _ in range(rng.choice([3, 4, 5]))] for _ in L2]
rk = [rank_mod(Y) for Y in famV]
dV = [[rk[a] + rk[b] - rank_mod(famV[a] + famV[b]) for b in range(len(L2))] for a in range(len(L2))]
dV = [[rk[a] - x for x in row] for a, row in enumerate(dV)]   # delta = dim Y_a - dim(Y_a cap Y_b)
for name, dl in (("coordinate", dC), ("random subspaces mod p", dV)):
    ab = {}
    for a in range(len(L2)):
        for b in range(len(L2)):
            s = dist(2, L2[a], L2[b]); t, c = ab.get(s, (0, 0)); ab[s] = (t + dl[a][b], c + 1)
    abar2 = {s: Fr(t, c) for s, (t, c) in ab.items()}
    up = {}
    for u in range(len(L3)):
        for v in range(len(L3)):
            r = dist(3, L3[u], L3[v])
            tot = sum(dl[proj[(u, p)]][proj[(v, p)]] for p in range(1, 64))
            t, c = up.get(r, (0, 0)); up[r] = (t + tot, c + 1)
    abar3 = {r: Fr(t, c) for r, (t, c) in up.items()}
    pred = {r: sum(c * abar2.get(s, 0) for s, c in closed(3, r).items()) for r in range(0, 4) if r > 0}
    ok = all(abar3[r] == pred[r] for r in (1, 2, 3))
    print(f"{name}: abar2 = {[str(abar2[s]) for s in (1, 2)]}; explicit abar3 = {[str(abar3[r]) for r in (1, 2, 3)]}; "
          f"recursion = {[str(pred[r]) for r in (1, 2, 3)]}; equal: {ok}")
    # is the explicit projected family distance-determined pairwise? (not needed; information only)
