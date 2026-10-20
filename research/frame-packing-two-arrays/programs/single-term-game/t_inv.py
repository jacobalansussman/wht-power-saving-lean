#!/usr/bin/env python3
"""t_inv.py: test candidate invariants for a proof of Conjecture T on every reachable state (exhaustive, stg.py):
 (C0) L1 + L2 <= M
 (C1) for every locked copy E = (x, y): L1 + lv(y) <= M and lv(x) + L2 <= M
 (C2) for locked copies E, F of DIFFERENT lineages: lv(x_E) + lv(y_F) <= M
 (C2s) same, E, F distinct copies of ONE lineage (expected to fail: branches share their prefix)
 (C3) for every locked copy: lv(x) + lv(y) <= M + (number of copies of its lineage) ... reported as max excess
Prints the maximum excess (left side minus M) of each over every state of each depth."""
import sys, time
from stg import *
from sweep import GAMES
from t_stg_toys import paid_classes

def excesses(g, st, M):
    A1, A2, lins = st
    L1, L2 = maxlv(g, A1), maxlv(g, A2)
    e = {"C0": L1 + L2 - M, "C1": -99, "C2": -99, "C2s": -99}
    copies = []
    for li, (t, xs) in enumerate(lins):
        for x in bits(xs):
            y = g.tw[t][x]; copies.append((li, g.lv[x], g.lv[y]))
    for (li, a, b) in copies:
        e["C1"] = max(e["C1"], L1 + b - M, a + L2 - M)
    for (li, a, b) in copies:
        for (lj, c, d) in copies:
            if (li, a, b) == (lj, c, d): continue
            k = "C2" if li != lj else "C2s"
            e[k] = max(e[k], a + d - M)
    # position-level Phi = sum_k rank of the constraints at level k (an UPPER bound on the module's
    # sum_k dim eps(G_k): every lineage entry treated as an independent generic ratio)
    N = g.N
    adjp = {}
    for li, (t, xs) in enumerate(lins):
        for x in bits(xs):
            u, w = x, N + g.tw[t][x]
            adjp.setdefault(u, []).append((w, li)); adjp.setdefault(w, []).append((u, li))
    comps = []; seenv = set()
    for st0 in adjp:
        if st0 in seenv: continue
        mono = {st0: (1, ())}; q = [st0]
        while q:
            v = q.pop()
            sg, mo = mono[v]
            for (nb, li) in adjp[v]:
                if nb in mono: continue
                d = dict(mo)
                # psi(2,y) = -psi(1,x)/r ; psi(1,x) = -r psi(2,y)
                d[li] = d.get(li, 0) + (-1 if v < N else 1)
                mono[nb] = (-sg, tuple(sorted((k2, v2) for k2, v2 in d.items() if v2)))
                q.append(nb)
        seenv |= set(mono); comps.append(mono)
    Phi = 0
    for k in range(1, g.m + 1):
        dirs = set()
        if any(g.lv[x] >= k for x in bits(A1)): dirs.add("e1")
        if any(g.lv[y] >= k for y in bits(A2)): dirs.add("e2")
        for mono in comps:
            H1 = {mono[v] for v in mono if v < N and g.lv[v] >= k}
            H2 = {mono[v] for v in mono if v >= N and g.lv[v - N] >= k}
            f1 = len(H1) >= 2; f2 = len(H2) >= 2
            if (f1 and H2) or (f2 and H1): dirs |= {"e1", "e2"}; continue
            if f1: dirs.add("e1")
            if f2: dirs.add("e2")
            if H1 and H2 and not f1 and not f2:
                (s1, m1), = H1; (s2, m2), = H2
                d = dict(m1)
                for k2, v2 in m2: d[k2] = d.get(k2, 0) - v2
                dirs.add(("line", s1 * s2, tuple(sorted((k2, v2) for k2, v2 in d.items() if v2))))
        Phi += min(2, len(dirs))
    e["Phi"] = Phi - M
    # (S): every input-1 vertex u and input-2 vertex w that are witnessed by two DIFFERENT objects (available single,
    # or copies of different lineages, or u, w in one locked component at path distance >= 3): lv(u) + lv(w) <= M
    own1 = {}; own2 = {}
    for x in bits(A1): own1.setdefault(x, set()).add("A")
    for y in bits(A2): own2.setdefault(y, set()).add("A")
    for li, (t, xs) in enumerate(lins):
        for x in bits(xs):
            own1.setdefault(x, set()).add(li); own2.setdefault(g.tw[t][x], set()).add(li)
    e["S"] = -99
    for u, o1 in own1.items():
        for w, o2 in own2.items():
            if len(o1 | o2) >= 2 or ("A" in o1 and "A" in o2):
                e["S"] = max(e["S"], g.lv[u] + g.lv[w] - M)
    # component-level tests
    N = g.N
    adj = {}
    for li, (t, xs) in enumerate(lins):
        for x in bits(xs):
            u, w = x, N + g.tw[t][x]
            adj.setdefault(u, []).append((w, li)); adj.setdefault(w, []).append((u, li))
    e["C1c"] = -99; e["C2c"] = -99; e["C2cs"] = -99
    seen = set()
    for start in adj:
        if start in seen: continue
        comp = {start}; st = [start]
        while st:
            v = st.pop()
            for (nb, li) in adj[v]:
                if nb not in comp: comp.add(nb); st.append(nb)
        seen |= comp
        a = max(g.lv[v] for v in comp if v < N); b = max(g.lv[v - N] for v in comp if v >= N)
        e["C1c"] = max(e["C1c"], L1 + b - M, a + L2 - M)
        # paths of length >= 3 from each side-1 vertex (tree): BFS with depth and first-edge lineage
        for u in [v for v in comp if v < N]:
            dist = {u: (0, None)}; q = [u]
            while q:
                v = q.pop(0)
                d, fl = dist[v]
                for (nb, li) in adj[v]:
                    if nb not in dist:
                        dist[nb] = (d + 1, li if fl is None else fl); q.append(nb)
                        if nb >= N and d + 1 >= 3:
                            # last edge lineage li, first edge lineage dist[nb][1]
                            k = "C2c" if dist[nb][1] != li else "C2cs"
                            e[k] = max(e[k], g.lv[u] + g.lv[nb - N] - M)
    return e

def run(g, D, cap=300000):
    st = closure(g, 1, 1, ())
    layer = {canon(g, *st): st}
    ex_states = {}
    for M in range(1, D + 1):
        new = {}
        for key, s in layer.items():
            for (kind, s2) in successors(g, s):
                s2 = closure(g, *s2); k2 = canon(g, *s2)
                if k2 not in new: new[k2] = s2
            if len(new) > cap: break
        layer = new
        worst = {"C0": -99, "C1": -99, "C2": -99, "C2s": -99, "C1c": -99, "C2c": -99, "C2cs": -99, "S": -99, "Phi": -99}
        for k, s in layer.items():
            e = excesses(g, s, M)
            for kk in worst:
                if e[kk] > worst[kk]:
                    worst[kk] = e[kk]
                    if e[kk] > 0 and kk not in ex_states: ex_states[kk] = (M, s)
        print("  %s M=%d states=%d  max excess: %s%s" % (g.name, M, len(layer), "  ".join("%s %+d" % (k, v) for k, v in worst.items()),
              "  (cap)" if len(layer) > cap else ""), flush=True)
        if len(layer) > cap: break
    for kk, (M, s) in ex_states.items():
        A1, A2, lins = s
        print("   first %s failure at M=%d: A1 levels %s, A2 levels %s, lineages %s" % (kk, M,
              sorted(g.lv[x] for x in bits(A1)), sorted(g.lv[x] for x in bits(A2)),
              [(t, [(g.lv[x], g.lv[g.tw[t][x]]) for x in bits(xs)]) for (t, xs) in lins]), flush=True)

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        nm, ci, D = spec.split(":")
        if nm.startswith("H"): g = hamming_game(int(nm[1]))
        else:
            G, P, _ = GAMES[nm](); g = from_toylab(G, P, paid_classes(G, P)[int(ci)], name=nm + "/c" + ci)
        print("==", g.name, "m", g.m, flush=True)
        run(g, int(D))
