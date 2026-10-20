#!/usr/bin/env python3
# frame_graphs.py : the frame graphs (dual polar graphs of Sp(2m,2)) at block sizes 2,3,4, built from scratch:
# vertices = m-dimensional subspaces of the 2m-bit vectors on which the pairing x.Jy vanishes; adjacent = share m-1 dimensions.
# Prints the numbers an embedding of a code graph would have to respect; and the engine's exponent for the translated unit.
import itertools, collections, math, sys
def build(m):
    n = 2 * m
    def form(x, y):                      # symplectic pairing: bit k of x with bit k+m of y, plus the reverse
        a = (x & ((1 << m) - 1)) & (y >> m); b = (x >> m) & (y & ((1 << m) - 1))
        return (bin(a).count("1") + bin(b).count("1")) & 1
    def span(gens):
        s = {0}
        for g in gens: s |= {x ^ g for x in s}
        return frozenset(s)
    base = span([1 << k for k in range(m)])
    seen = {base: 0}; order = [base]; adj = collections.defaultdict(set); fr = [base]
    while fr:
        nx = []
        for U in fr:
            Ul = sorted(U)
            # hyperplanes of U = kernels of nonzero functionals; enumerate via bases
            hyps = set()
            basis = []
            sp = {0}
            for x in Ul:
                if x not in sp: basis.append(x); sp |= {y ^ x for y in sp}
            for mask in range(1, 1 << m):  # functional takes value bit k of mask on basis[k]
                H = frozenset(x for x in U if sum(((mask >> k) & 1) for k in range(m) if _coef(x, basis, k)) % 2 == 0)
                hyps.add(H)
            for H in hyps:
                perp = [v for v in range(1 << n) if all(form(v, h) == 0 for h in H)]
                for v in perp:
                    if v in U: continue
                    V = frozenset(H | {h ^ v for h in H})
                    if V not in seen: seen[V] = len(order); order.append(V); nx.append(V)
                    adj[seen[U]].add(seen[V])
        fr = nx
    return order, adj
_cache = {}
def _coef(x, basis, k):
    # coefficient of basis[k] in x (x in span(basis)); brute force over subsets, cached
    key = (x, tuple(basis))
    if key not in _cache:
        for mask in range(1 << len(basis)):
            s = 0
            for t in range(len(basis)):
                if (mask >> t) & 1: s ^= basis[t]
            _cache[(s, tuple(basis))] = mask
    return (_cache[key] >> k) & 1
for m in (2, 3, 4):
    order, adj = build(m)
    N = len(order); degs = collections.Counter(len(adj[v]) for v in range(N))
    d = {0: 0}; fr = [0]
    while fr:
        nx = []
        for u in fr:
            for w in adj[u]:
                if w not in d: d[w] = d[u] + 1; nx.append(w)
        fr = nx
    tri = sum(1 for u in range(N) for v in adj[u] if v > u for w in adj[v] if w > v and w in adj[u])
    cn = collections.Counter()
    for v in range(1, N):
        if d[v] <= 2: cn[(d[v], len(adj[0] & adj[v]))] += 1
    print(f"block size {m}: frames {N}, degree {sorted(degs.items())}, frames at distance k from one frame {sorted(collections.Counter(d.values()).items())}, "
          f"triangles {tri}, edges {sum(len(adj[v]) for v in range(N)) // 2}, (distance, common neighbours): count {sorted(cn.items())}, bipartite: {tri == 0}")
    sys.stdout.flush()
print("block size 5 by the same formulas [derived, not built]: frames 3*5*9*17*33 =", 3 * 5 * 9 * 17 * 33, ", degree 62, opposite frames per frame 2^15 =", 2 ** 15)
print("--- exponent of the translated unit by the engine's inequality in the limit of large s: sum over blocks (k/m)^z < W")
def root(ranks_count, m, W):             # ranks_count: {rank: number of blocks}
    f = lambda z: sum(c * (k / m) ** z for k, c in ranks_count.items()) - W
    if f(1.0) >= 0: return None
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return hi
for name, W, m, moves in [("Braverman-He causal code as it stands", 157, 5, 735), ("their nondeterministic certificate (not a schedule)", 166, 5, 744)]:
    z = root({1: moves}, m, W); rho = 1 - moves / (W * m)
    print(f"{name}: W={W}, m={m}, moves={moves} against {W*m}; saved fraction {rho:.4f}; all blocks of rank 1: x = {1-z:.4f};"
          f" per-rank formula -ln(1-rho)/ln m = {-math.log(1-rho)/math.log(m):.4f}; whole-block bound rho/((1-rho) ln m) = {rho/((1-rho)*math.log(m)):.4f};"
          f" same bound if every move sat in a rank-(m-1) block (not available to this code): {rho/((1-rho)*math.log(m/(m-1))):.4f}")
for m in (3, 4, 5, 7):
    print(f"saved fraction needed for x = 0.1 with rank-1 blocks at block size {m}: {1 - m ** -0.1:.4f}  (moves per array at most {m * m ** -0.1:.4f} of {m})")
print("rate a code of this design needs (K sessions on a graph with E edges before deletion, ends of a session m apart): cost (E-K)/(mK);")
for m in (3, 5):
    print(f"   m={m}: gain at all needs K/E > {1/(m+1):.4f}; x >= 0.1 needs K/E > {1/(1 + m * m ** -0.1):.4f};  the preprint's figures: 166/910 = {166/910:.4f} (certificate), 157/(735+157) = {157/892:.4f} (causal)")
