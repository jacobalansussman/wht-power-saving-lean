# quot3.py -- ceiling curve rho_max(k) for contractions of Gamma_m, by two constructive searches:
#  (A) agglomerative: start from singletons (rho = m), repeatedly merge the adjacent pair of parts that lowers rho
#      least; record the best rho found at every k.
#  (B) distance-profile partitions: class of F = (d(F, r_1), ..., d(F, r_j)) for random reference frames r_i.
# Both give ACHIEVABLE rho (lower estimates of rho_max(k)); together with quot2's annealing they map how coarse a
# contraction can be while still certifying a given number of moves per session.
import sys, random, json, time, collections
import numpy as np
import dpg

def setup(m):
    verts, idx, adj, dist = dpg.build(m)
    n = len(verts)
    Dm = np.array(dist, dtype=np.int64)
    OPP = (Dm == m).astype(np.int64)
    ADJ = (Dm == 1).astype(np.int64)
    return n, Dm, OPP, ADJ

def rho_of(P, C, D, npairs, nedges):
    k = P.shape[0]
    off = ~np.eye(k, dtype=bool)
    cross = C[off].sum() / 2
    return ((P * D).sum() / npairs) / (cross / nedges)

def qmetric(C):
    k = C.shape[0]
    D = np.where(C > 0, 1, 10 ** 6).astype(np.int64); np.fill_diagonal(D, 0)
    for w in range(k): D = np.minimum(D, D[:, [w]] + D[[w], :])
    return D

def agglomerate(m, n, Dm, OPP, ADJ, record):
    npairs = OPP.sum(); nedges = ADJ.sum() // 2
    P = OPP.copy(); C = ADJ.copy(); D = Dm.copy()
    members = [[i] for i in range(n)]
    curve = {n: float(rho_of(P, C, D, npairs, nedges))}
    t0 = time.time()
    while P.shape[0] > 2:
        k = P.shape[0]
        off = ~np.eye(k, dtype=bool)
        cross = C[off].sum() / 2
        num = (P * D).sum()
        best = None
        ai, bi = np.nonzero(np.triu(C, 1))
        for a, b in zip(ai, bi):
            # merged metric (before deleting row b): D' = min(D, D[:,a]+D[b,:], D[:,b]+D[a,:])
            Dn = np.minimum(D, np.minimum(D[:, [a]] + D[[b], :], D[:, [b]] + D[[a], :]))
            # pairs inside the merged part get distance 0: rows/cols a,b
            numn = (P * Dn).sum() - 2 * P[a, b] * Dn[a, b]
            # after merge, (a,b) pairs are internal: Dn[a,b] should be 0 -- min formula gives D[a,a]+D[b,b]... = 0 already
            crossn = cross - C[a, b]
            if crossn <= 0: continue
            val = (numn / npairs) / (crossn / nedges)
            if best is None or val > best[0]: best = (val, a, b)
        val, a, b = best
        # perform merge of b into a
        Dn = np.minimum(D, np.minimum(D[:, [a]] + D[[b], :], D[:, [b]] + D[[a], :]))
        P[a, :] += P[b, :]; P[:, a] += P[:, b]
        C[a, :] += C[b, :]; C[:, a] += C[:, b]
        keep = [i for i in range(k) if i != b]
        P = P[np.ix_(keep, keep)]; C = C[np.ix_(keep, keep)]; D = Dn[np.ix_(keep, keep)]
        members[a] = members[a] + members[b]; del members[b]
        # recompute exactly (guards against any slip in the incremental formula)
        Dq = qmetric(np.where(np.eye(len(keep), dtype=bool), 0, C))
        assert (Dq == D).all(), "metric update mismatch"
        r = float(rho_of(P, C, D, npairs, nedges))
        curve[len(keep)] = r
        if len(keep) in record:
            print(f"   agglomerative m={m} k={len(keep)}: rho = {r:.5f}  sizes {sorted(len(x) for x in members)[-6:]} ({time.time()-t0:.0f}s)")
            sys.stdout.flush()
    return curve

def profiles(m, n, Dm, OPP, ADJ, tries, seed, kmax=10 ** 9):
    rng = random.Random(seed)
    OPPf = OPP.astype(float); ADJf = ADJ.astype(float)
    npairs = OPP.sum(); nedges = ADJ.sum() // 2
    best = {}
    for t in range(tries):
        j = rng.randint(1, 6)
        R = [rng.randrange(n) for _ in range(j)]
        keys = [tuple(int(Dm[F, r]) for r in R) for F in range(n)]
        cls = {}; part = np.array([cls.setdefault(kk, len(cls)) for kk in keys])
        k = len(cls)
        if k > kmax: continue
        M = np.zeros((k, n)); M[part, np.arange(n)] = 1.0
        P = np.rint(M @ OPPf @ M.T).astype(np.int64); C = np.rint(M @ ADJf @ M.T).astype(np.int64)
        D = qmetric(np.where(np.eye(k, dtype=bool), 0, C))
        r = float(rho_of(P, C, D, npairs, nedges))
        if k not in best or r > best[k][0]: best[k] = (r, j, [int(x) for x in sorted(collections.Counter(part.tolist()).values())])
    return best

if __name__ == '__main__':
    m = int(sys.argv[1]); mode = sys.argv[2]
    n, Dm, OPP, ADJ = setup(m)
    print(f"m={m}: {n} frames, {ADJ.sum()//2} edges, {OPP.sum()} ordered opposite pairs"); sys.stdout.flush()
    out = {}
    if 'prof' in mode:
        best = profiles(m, n, Dm, OPP, ADJ, int(sys.argv[3]), 1)
        print("  distance-profile partitions, best rho per k:")
        for k in sorted(best): print(f"   k={k}: rho = {best[k][0]:.5f} (refs {best[k][1]}, sizes {best[k][2][-8:]})")
        out['profiles'] = {k: best[k][0] for k in best}
    if 'agg' in mode:
        rec = set([2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 134])
        curve = agglomerate(m, n, Dm, OPP, ADJ, rec)
        out['agglomerative'] = curve
    json.dump(out, open(f"quot3_m{m}_{mode}.json", "w"))
