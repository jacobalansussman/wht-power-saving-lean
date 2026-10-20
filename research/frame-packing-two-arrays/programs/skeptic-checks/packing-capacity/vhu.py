#!/usr/bin/env python3
"""vhu.py (v-fp-capacity): Theorem A(b) chain on rebuilt W=2 units, own code: in the raw gauge and random gauges,
capacities c(e) = number of live virtual moves (Sp-label classes) along each edge of Gamma_m; check that the
two-commodity fractional flow S_a -> T_a (value 1 each) is feasible (Hu), that every virtual move is an edge, and
that sum_j k_j >= 2m. Negative control: capacity 1 per paid move placed on ONE live class's edge (the first)."""
import sys, json, random
import numpy as np
from scipy.optimize import linprog
from vengine import Eng, run, analyse
from vunits import words

def lag_dist(m, F1, F2): return m - (len(F1 & F2).bit_length() - 1)

def hu_feasible(edges_cap, S, T):
    nodes = sorted({x for e in edges_cap for x in e} | set(S) | set(T), key=lambda z: sorted(z))
    idx = {x: i for i, x in enumerate(nodes)}
    arcs = []
    for (u, v) in edges_cap: arcs += [(idx[u], idx[v]), (idx[v], idx[u])]
    na = len(arcs); nv = 2 * na
    Aeq = []; beq = []
    for a in range(2):
        for x in range(len(nodes)):
            row = np.zeros(nv)
            for i, (u, v) in enumerate(arcs):
                if u == x: row[a * na + i] += 1
                if v == x: row[a * na + i] -= 1
            rhs = (1 if x == idx[S[a]] else 0) - (1 if x == idx[T[a]] else 0)
            Aeq.append(row); beq.append(rhs)
    Aub = []; bub = []
    for k, (e, c) in enumerate(edges_cap.items()):
        row = np.zeros(nv)
        for a in range(2): row[a * na + 2 * k] = 1; row[a * na + 2 * k + 1] = 1
        Aub.append(row); bub.append(c)
    if not Aub: return S == T
    r = linprog(np.zeros(nv), A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=(0, None), method="highs")
    return r.status == 0

def main(seed, nr):
    rng = random.Random(seed); rows = []
    for m in (2, 3, 4):
        E = Eng(m); L0 = E.L0()
        letters = [E.dirm(z) for z in range(1, E.n)] + [E.phase(z) for z in range(1, E.n)] + [E.Xm(z) for z in range(1, E.n)]
        Kinv = E.spinv(E.Kk)
        ident = tuple(1 << j for j in range(2 * m))
        for (name, W, R, w, dl) in words(E, rng, nr):
            if W != 2: continue
            snaps = run(E, W, R, w)
            gauges = [("raw", [ident, ident], [ident, ident])]
            for gi in range(3):
                hs = []; hinv = []
                for a in range(2):
                    M = E.I
                    for _ in range(10): M = letters[rng.randrange(len(letters))] @ M
                    M = np.round(M * (1 << 24)) / (1 << 24)
                    hs.append(E.sp_of_mat(M)); hinv.append(E.sp_of_mat(M.conj().T))
                gauges.append(("random%d" % gi, hs, hinv))
            for (gname, hs, hinv) in gauges:
                cls = lambda g, l, hs=hs: E.comp(E.sp(g), hs[l])
                nodes, live = analyse(E, W, R, w, snaps, dl, cls)
                pos = lambda g, l, hinv=hinv: E.frame(E.comp(hinv[l], E.spinv(g)), L0)
                cap = {}; cap1 = {}; nonedge = 0; sumk = 0; M_ = 0
                for t, st in enumerate(w):
                    if st[0] != "move": continue
                    M_ += 1; r = st[1]; seen = {}
                    for (g, l) in snaps[t][r]:
                        k = cls(g, l)
                        if k in live[t][r] and k not in seen:
                            g2 = E.key(st[2] @ E.mats[g]); seen[k] = (pos(g, l), pos(g2, l))
                    sumk += len(seen); first = True
                    for k, (X, Y) in seen.items():
                        if lag_dist(m, X, Y) != 1: nonedge += 1
                        e = tuple(sorted((X, Y), key=lambda z: sorted(z)))
                        cap[e] = cap.get(e, 0) + 1
                        if first: cap1[e] = cap1.get(e, 0) + 1; first = False
                S = [E.frame(hinv[a], L0) for a in range(2)]
                T = [E.frame(E.comp(hinv[a], Kinv), L0) for a in range(2)]
                assert all(lag_dist(m, S[a], T[a]) == m for a in range(2))
                rec = dict(m=m, name=name, gauge=gname, M=M_, sum_k=sumk, nonedge=nonedge,
                           hu=hu_feasible(cap, S, T), hu_one_per_move=hu_feasible(cap1, S, T))
                rec["ok"] = rec["hu"] and nonedge == 0 and sumk >= 2 * m
                rows.append(rec); print(json.dumps(rec), flush=True)
    bad = [r for r in rows if not r["ok"]]; neg = [r for r in rows if not r["hu_one_per_move"]]
    print("SUMMARY rows=%d failures=%d negative-control-infeasible=%d" % (len(rows), len(bad), len(neg)))
    json.dump(rows, open("vhu_s%d.json" % seed, "w"), indent=1)

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
