#!/usr/bin/env python3
"""cx-medium / cliquefloor.py (own code, stdlib; python3 -I).  A lower bound on the centre charge of ONE call, valid for
every decoder (any real gather/scatter coefficients, any number of totals), in the hub's accounting (gcert E2: the
targets wait at frame 0 at the scatter; a total's copy stands at a frame containing every gathered label, so it costs
at least c(S) = dim span_F2{u_b : b in S}, S its gather set).
Argument [proved here, short]: let C be a clique of omega ports with labels pairwise NON-orthogonal (B = 1 on C x C,
i.e. mutually orthogonal roots).  The delivered map must be the identity on C x C (no side information can pass inside
C), and each total contributes a rank <= 1 piece there, nonzero only if its gather set S meets C (and its scatter set
does).  So at least omega totals have gather sets meeting C.  Summing over all maximum cliques:
    sum_k T(S_k) >= omega * #cliques,   T(S) = number of maximum cliques that S meets.
Each total costs c(S_k) >= c(W_k) with W_k = span of its labels, and T(S_k) <= T(ports with label in W_k).  Hence
    centre charge per call >= omega * #cliques * min_W dim(W) / T(labels in W)
over all subspaces W of the label space F_2^k (enumerated exhaustively here: 29,212 for k = 7, 374 for k = 5)."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from designs_cx import design, labels, dot

def subspaces(k):
    """all subspaces of F_2^k as frozensets of elements, by closure from spans (BFS on dimension)"""
    seen = {frozenset([0])}
    layer = [frozenset([0])]
    out = [frozenset([0])]
    for d in range(k):
        nxt = []
        for S in layer:
            for v in range(1, 1 << k):
                if v in S: continue
                T = frozenset(S | {x ^ v for x in S})
                if T not in seen:
                    seen.add(T); nxt.append(T); out.append(T)
        layer = nxt
    return out

def max_cliques(n, adj, size):
    res = []
    def ext(R, P):
        if len(R) == size: res.append(tuple(R)); return
        for v in sorted(P):
            if v > (R[-1] if R else -1):
                ext(R + [v], P & adj[v])
    ext([], set(range(n)))
    return res

out = {}
for name in sys.argv[1:] or ['E7', 'D6', 'A5']:
    L = design(name); labs, k, B = labels(L); n = len(L)
    adj = [set(j for j in range(n) if j != i and B[i][j] == 1) for i in range(n)]
    # omega: the largest clique size (brute force upward)
    om = 1
    while max_cliques(n, adj, om + 1): om += 1
    cl = max_cliques(n, adj, om)
    of = {}
    for ci, c in enumerate(cl):
        for p in c: of.setdefault(p, set()).add(ci)
    lab2ports = {}
    for p, u in enumerate(labs): lab2ports.setdefault(u, []).append(p)
    best = None
    for W in subspaces(k):
        d = len(W).bit_length() - 1
        if d == 0: continue
        ports = [p for u in W if u in lab2ports for p in lab2ports[u]]
        T = len(set().union(*[of.get(p, set()) for p in ports])) if ports else 0
        if T == 0: continue
        r = d / T
        if best is None or r < best[0]: best = (r, d, T, len(ports))
    need = om * len(cl)
    floor = need * best[0]
    out[name] = dict(n=n, k=k, omega=om, cliques=len(cl), incidences_needed=need, best_ratio=best[0], best_dim=best[1],
                     best_T=best[2], best_ports=best[3], floor_per_call=floor, floor_per_port=floor / n,
                     break_even_P1=2 * n / 3, break_even_P2=4 * n / 5)
    print(f"{name}: n={n} k={k} omega={om} max cliques={len(cl)} -> totals must meet cliques >= {need} times; "
          f"best dim(W)/T = {best[1]}/{best[2]} (W holds {best[3]} ports) -> centre >= {floor:.2f} per call "
          f"({floor/n:.3f} per port); break-even 1 pair {2*n/3:.1f}, 2 pairs {4*n/5:.1f}", flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cliquefloor.json'), 'w'), indent=1)
