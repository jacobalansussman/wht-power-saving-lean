# quot2.py -- incremental annealing of the best k-part contraction of Gamma_m (m = 2, 3, 4).
# rho(pi) = E_{ordered opposite pairs}[d_q] / (share of edges crossing pi): the ROUTING value of the quotient network
# for the uniform opposite-pair demand, per session. Any lower bound computed on the quotient (cuts, Shannon,
# Ingleton, DFZ, non-Shannon: anything valid for codes on the quotient) is <= rho(pi), because routing is a code.
# So max_pi rho(pi) over k-part pi is a CEILING on every quotient method with k parts.
import sys, random, math, json, time, collections
import numpy as np
import dpg

def qmetric(C, k):
    INF = 10 ** 6
    D = np.full((k, k), INF, dtype=np.int64)
    np.fill_diagonal(D, 0)
    D[(C > 0)] = 1
    np.fill_diagonal(D, 0)
    for w in range(k):
        D = np.minimum(D, D[:, [w]] + D[[w], :])
    return D

class State:
    def __init__(self, n, k, adj, opp, part, npairs, nedges):
        self.n, self.k, self.adj, self.opp = n, k, adj, opp
        self.part = list(part); self.npairs = npairs; self.nedges = nedges
        self.P = np.zeros((k, k), dtype=np.int64); self.C = np.zeros((k, k), dtype=np.int64)
        for s in range(n):
            for t in opp[s]: self.P[part[s], part[t]] += 1
            for t in adj[s]: self.C[part[s], part[t]] += 1   # each edge counted twice (both ends)
    def value(self):
        k = self.k
        Cx = self.C.copy(); np.fill_diagonal(Cx, 0)
        cross = Cx.sum() / 2
        if cross == 0: return 0.0
        D = qmetric(Cx, k)
        if (D >= 10 ** 6).any() and (self.P[D >= 10 ** 6] > 0).any(): return 0.0
        num = (self.P * np.where(D >= 10 ** 6, 0, D)).sum()
        return (num / self.npairs) / (cross / self.nedges)
    def move(self, v, b):
        a = self.part[v]
        if a == b: return
        for t in self.opp[v]:
            pt = self.part[t]
            self.P[a, pt] -= 1; self.P[pt, a] -= 1; self.P[b, pt] += 1; self.P[pt, b] += 1
        for t in self.adj[v]:
            pt = self.part[t]
            self.C[a, pt] -= 1; self.C[pt, a] -= 1; self.C[b, pt] += 1; self.C[pt, b] += 1
        self.part[v] = b

def anneal(G, m, k, iters, seed, T0=0.05):
    verts, idx, adj, dist = G
    n = len(verts)
    opp = [[t for t in range(n) if dist[s][t] == m] for s in range(n)]
    npairs = sum(len(o) for o in opp); nedges = sum(len(a) for a in adj) // 2
    rng = random.Random(seed)
    st = State(n, k, adj, opp, [rng.randrange(k) for _ in range(n)], npairs, nedges)
    cur = st.value(); best = (cur, st.part[:])
    T = T0
    for it in range(iters):
        v = rng.randrange(n); a = st.part[v]
        # propose: move to the part of a random neighbour (keeps parts connected-ish) or random
        b = st.part[rng.choice(adj[v])] if rng.random() < 0.7 else rng.randrange(k)
        if b == a: continue
        st.move(v, b)
        val = st.value()
        if val >= cur or rng.random() < math.exp((val - cur) / T):
            cur = val
            if val > best[0] + 1e-12: best = (val, st.part[:])
        else:
            st.move(v, a)
        T = max(1e-4, T * (1 - 5.0 / iters))
    return best, opp, npairs, nedges

if __name__ == '__main__':
    m = int(sys.argv[1]); ks = [int(x) for x in sys.argv[2].split(',')]; iters = int(sys.argv[3]); restarts = int(sys.argv[4])
    t0 = time.time()
    G = dpg.build(m)
    print(f"m={m}: built {len(G[0])} frames in {time.time()-t0:.1f}s"); sys.stdout.flush()
    out = {}
    for k in ks:
        res = []
        for r in range(restarts):
            best, opp, npairs, nedges = anneal(G, m, k, iters, 7919 * k + r)
            res.append(best)
        res.sort(key=lambda t: -t[0])
        val, part = res[0]
        sizes = sorted(collections.Counter(part).values())
        print(f"m={m} k={k}: best rho = {val:.5f}  top: {[round(x[0], 4) for x in res[:4]]}  parts used {len(sizes)} sizes {sizes}  ({time.time()-t0:.0f}s)")
        sys.stdout.flush()
        out[k] = dict(rho=val, part=part, top=[x[0] for x in res])
    json.dump(out, open(f"quot2_m{m}.json", "w"))
