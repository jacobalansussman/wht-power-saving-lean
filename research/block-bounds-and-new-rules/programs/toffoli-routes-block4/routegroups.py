"""routegroups.py m -- every shortest route to the kernel as a factorization K = Q (a+b t_m)...(a+b t_1) h_0
(t_j involutive monomial matrices, Q and h_0 monomial). For every route through a representative of each AGL-orbit
of first backward steps (census_m{m}.json), compute the involutions t_j in one common frame and the group they
generate (permutation part and with phases). Invariants are unchanged by the route symmetry Lambda (right
multiplication), so weighting by orbit sizes gives the distribution over all routes."""
import sys
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, time, collections
import numpy as np
from rt import Game, A_, B_

m = int(sys.argv[1])
g = Game(m); n = g.n
cen = json.load(open("census_m%d.json" % m)  )

def mono(M):
    """monomial matrix -> (perm, exps) with M e_x = i^exps[x] e_perm[x]"""
    A = np.abs(M) > 1e-9
    assert (A.sum(axis=0) == 1).all() and (A.sum(axis=1) == 1).all(), "not monomial"
    p = A.argmax(axis=0)
    v = M[p, np.arange(n)]
    e = np.zeros(n, dtype=np.int64)
    for k, z in enumerate([1, 1j, -1, -1j]):
        e[np.abs(v - z) < 1e-9] = k
    assert all(min(abs(x - z) for z in [1, 1j, -1, -1j]) < 1e-9 for x in v), "phase not in Z4"
    return tuple(p.tolist()), tuple(e.tolist())

def mat(t):
    p, e = t
    M = np.zeros((n, n), dtype=complex)
    for x in range(n):
        M[p[x], x] = 1j ** e[x]
    return M

def comp(t, r):  # t o r
    pt, et = t; pr, er = r
    return (tuple(pt[pr[x]] for x in range(n)), tuple((er[x] + et[pr[x]]) % 4 for x in range(n)))

def closure(gens, cap=200000):
    idt = (tuple(range(n)), tuple([0] * n))
    seen = {idt}; frontier = [idt]
    while frontier:
        nf = []
        for x in frontier:
            for s in gens:
                y = comp(s, x)
                if y not in seen:
                    seen.add(y); nf.append(y)
                    if len(seen) > cap:
                        return None
        frontier = nf
    return seen

def group_invariants(G):
    els = list(G)
    order = len(els)
    idt = (tuple(range(n)), tuple([0] * n))
    inv = sum(1 for x in els if x != idt and comp(x, x) == idt)
    ab = all(comp(x, y) == comp(y, x) for x in els[:64] for y in els[:64]) and \
        all(comp(x, y) == comp(y, x) for x in els for y in els) if order <= 64 else None
    center = sum(1 for x in els if all(comp(x, y) == comp(y, x) for y in els)) if order <= 64 else None
    expo = 1
    for x in els:
        k = 1; y = x
        while y != idt:
            y = comp(x, y); k += 1
        expo = max(expo, k)
    return dict(order=order, involutions=inv, abelian=ab, center=center, exponent=expo)

def routes_from(h, lev):
    """yield lists [h_lev, ..., h_0] of elements along every halving chain down to B"""
    if lev == 0:
        yield [h]; return
    for kk, (M, c) in g.children(h).items():
        if g.continues(c, lev - 1):
            for tail in routes_from(c, lev - 1):
                yield [h] + tail

def taus_of(chain):
    """chain = [h_m, ..., h_0]; returns t_1..t_m in h_0's row frame, checks the factorization"""
    hs = chain[::-1]   # h_0 .. h_m
    sig = np.zeros((n, n), dtype=complex)
    for k in range(n // 2):
        sig[2 * k + 1, 2 * k] = 1; sig[2 * k, 2 * k + 1] = 1
    Q = np.eye(n, dtype=complex)
    taus = []
    for j in range(1, m + 1):
        fwd = g.snap((A_ * np.eye(n) + B_ * sig) @ hs[j - 1])
        P = g.snap(fwd @ np.conj(hs[j]).T)
        mono(P)
        tau = g.snap(Q @ sig @ np.conj(Q).T)
        taus.append(tau)
        Q = g.snap(Q @ P)
    # check: (a+b t_m)...(a+b t_1) h_0 == Q h_m
    prod = hs[0]
    for tau in taus:
        prod = g.snap((A_ * np.eye(n) + B_ * tau) @ prod)
    assert np.allclose(prod, Q @ hs[m]), "factorization check failed"
    return [mono(t) for t in taus]

t0 = time.time()
dist = collections.Counter()
dist_w = collections.Counter()
examples = {}
for r in cen:
    if not r["cont"]:
        continue
    M = [tuple(x) for x in r["rep"]]
    E0 = g.compat(g.K)
    h3 = g.back_step(g.K, M, E0)
    cnt = 0
    allr = list(routes_from(h3, m - 1))
    assert len(allr) == r["nroutes"]
    import random as _rnd
    _rnd.seed(7)
    SAMPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    pick = allr if len(allr) <= SAMPLE else _rnd.sample(allr, SAMPLE)
    wfac = r["size"] * len(allr) / len(pick)
    for tail in pick:
        chain = [g.K] + tail
        taus = taus_of(chain)
        perms = [(t[0], tuple([0] * n)) for t in taus]
        Gp = closure(perms, 70000)
        Gm = closure(taus, 5000)
        # sharply transitive: products t_S (descending order) send point 0 to distinct points
        imgs = set()
        for S in range(2 ** m):
            x = 0
            for j in range(m):
                if (S >> j) & 1:
                    x = taus[j][0][x]
            imgs.add(x)
        sharp = len(imgs) == n
        commute = all(comp(perms[i], perms[j]) == comp(perms[j], perms[i]) for i in range(m) for j in range(m))
        invp = group_invariants(Gp) if Gp else dict(order='>70000', abelian=None, exponent=None, involutions=None, center=None)
        key = (r["translation"], invp["order"], invp["abelian"], invp["exponent"], invp["involutions"],
               invp["center"], len(Gm) if Gm else None, sharp, commute)
        dist[key] += 1
        dist_w[key] += wfac
        examples.setdefault(key, [[list(t[0]) for t in taus], [list(t[1]) for t in taus]])
        cnt += 1
    print("rep size %d translation=%s routes %d, analysed %d t=%.1f" % (r["size"], r["translation"], len(allr), cnt, time.time() - t0), flush=True)
    assert cnt == len(pick)
print("key = (first step translation, |perm group|, abelian, exponent, #involutions, |center|, |monomial group|, sharply transitive, t_j commute)")
for k in sorted(dist, key=lambda k: str(k)):
    print(k, "routes from reps:", dist[k], " all routes (weighted):", dist_w[k])
print("total weighted", sum(dist_w.values()))
json.dump({str(k): v for k, v in examples.items()}, open("routegroups_examples_m%d.json" % m, "w"))
