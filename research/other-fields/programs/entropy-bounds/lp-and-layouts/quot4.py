# quot4.py -- refine the agglomerative contractions of quot3 by annealing (vertex moves between adjacent parts),
# to tighten the achievable rho at fixed k; also random-matching merges at block 4.
import sys, random, math, json, time, collections
import numpy as np
import dpg, quot2, quot3

def agglom_parts(m, n, Dm, OPP, ADJ, targets):
    """re-run agglomeration, returning the membership at each target k"""
    npairs = OPP.sum(); nedges = ADJ.sum() // 2
    P = OPP.copy(); C = ADJ.copy(); D = Dm.copy(); members = [[i] for i in range(n)]
    got = {}
    while P.shape[0] > 2:
        k = P.shape[0]; off = ~np.eye(k, dtype=bool); cross = C[off].sum() / 2
        best = None
        ai, bi = np.nonzero(np.triu(C, 1))
        for a, b in zip(ai, bi):
            Dn = np.minimum(D, np.minimum(D[:, [a]] + D[[b], :], D[:, [b]] + D[[a], :]))
            numn = (P * Dn).sum(); crossn = cross - C[a, b]
            if crossn <= 0: continue
            val = (numn / npairs) / (crossn / nedges)
            if best is None or val > best[0]: best = (val, a, b)
        val, a, b = best
        Dn = np.minimum(D, np.minimum(D[:, [a]] + D[[b], :], D[:, [b]] + D[[a], :]))
        P[a, :] += P[b, :]; P[:, a] += P[:, b]; C[a, :] += C[b, :]; C[:, a] += C[:, b]
        keep = [i for i in range(k) if i != b]
        P = P[np.ix_(keep, keep)]; C = C[np.ix_(keep, keep)]; D = Dn[np.ix_(keep, keep)]
        members[a] = members[a] + members[b]; del members[b]
        if len(keep) in targets:
            part = [0] * n
            for c, mem in enumerate(members):
                for v in mem: part[v] = c
            got[len(keep)] = part
    return got

def anneal_from(G, m, part, iters, seed, T0=0.01):
    verts, idx, adj, dist = G
    n = len(verts); k = max(part) + 1
    opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
    npairs = sum(len(o) for o in opp); nedges = sum(len(a) for a in adj) // 2
    rng = random.Random(seed)
    st = quot2.State(n, k, adj, opp, part, npairs, nedges)
    cur = st.value(); best = (cur, st.part[:]); T = T0
    for it in range(iters):
        v = rng.randrange(n); a = st.part[v]; b = st.part[rng.choice(adj[v])]
        if b == a: continue
        st.move(v, b)
        if st.part.count(a) == 0 and False: pass
        val = st.value()
        if val >= cur or rng.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0] + 1e-12: best = (val, st.part[:])
        else: st.move(v, a)
        T = max(1e-5, T * (1 - 5.0 / iters))
    return best

if __name__ == '__main__':
    m = int(sys.argv[1]); targets = [int(x) for x in sys.argv[2].split(',')]; iters = int(sys.argv[3])
    G = dpg.build(m)
    n, Dm, OPP, ADJ = quot3.setup(m)
    got = agglom_parts(m, n, Dm, OPP, ADJ, set(targets))
    out = {}
    for k in sorted(got, reverse=True):
        part = got[k]
        st0 = quot2.State(n, k, G[2], [[t for t in range(n) if G[3][s][t] == m] for s in range(n)], part,
                          int(OPP.sum()), int(ADJ.sum() // 2))
        v0 = st0.value()
        best = max((anneal_from(G, m, part, iters, s) for s in range(3)), key=lambda t: t[0])
        k_used = len(set(best[1]))
        print(f"m={m} k={k}: agglomerative rho {v0:.5f} -> annealed {best[0]:.5f} (parts used {k_used})"); sys.stdout.flush()
        out[k] = dict(agg=v0, annealed=best[0], parts_used=k_used)
    json.dump(out, open(f"quot4_m{m}.json", "w"))
