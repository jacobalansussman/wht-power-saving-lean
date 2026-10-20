"""clique_test.py: the tight-pair clique test (PROOF.md Theorem 1) for any commutative-slice target (m, S, W).

Tight forms: odd-trace binary quadratic forms Q with exactly W moves in their support.  Two moves are joined when
some tight form contains both; a move is 'live' when some tight form contains it.  Theorem 1: every clique of live
moves has at most W elements (K = B^T A has rank <= W and is unit triangular on such a clique).
Exact maximum clique by simple branch and bound (graphs here are small)."""
import sys
import numpy as np
from refl import quad_forms


def max_clique(adj, nodes):
    best = []
    def expand(R, P):
        nonlocal best
        if len(R) + len(P) <= len(best):
            return
        if not P:
            if len(R) > len(best):
                best = list(R)
            return
        for v in list(P):
            if len(R) + len(P) <= len(best):
                return
            expand(R + [v], [u for u in P if adj[v][u]])
            P = [u for u in P if u != v]
    expand([], list(nodes))
    return best


def clique_test(m, S, W, verbose=True, label=""):
    vals, tr = quad_forms(m)
    sub = vals[:, S]
    wt = sub.sum(axis=1)
    odd = tr == 1
    dmin = int(wt[odd].min())
    if dmin < W:
        if verbose:
            print(f"{label} m={m} W={W} M={len(S)}: deletion test FAILS (odd form of weight {dmin} < W)")
        return ("deletion", dmin)
    tight = np.nonzero(odd & (wt == W))[0]
    M = len(S)
    adj = [[False] * M for _ in range(M)]
    live = [False] * M
    for f in tight:
        D = np.nonzero(sub[f])[0]
        for a in D:
            live[a] = True
            for b in D:
                if a != b:
                    adj[a][b] = True
    nodes = [i for i in range(M) if live[i]]
    cl = max_clique(adj, nodes)
    verdict = "EXCLUDED" if len(cl) > W else "passes"
    if verbose:
        print(f"{label} m={m} W={W} M={len(S)}: tight odd forms {len(tight)}, live moves {len(nodes)}, "
              f"max clique {len(cl)} {[S[i] for i in cl]} -> {verdict}", flush=True)
    return (verdict, len(cl))


if __name__ == "__main__":
    cases = [
        ("plain 1in3 b3", 3, [1, 2, 4], 1),
        ("plain 2in6 b3", 3, [1, 1, 2, 2, 4, 4], 2),
        ("plain 3in9 b3", 3, [1, 2, 4] * 3, 3),
        ("plain 2in8 b4", 4, [1, 2, 4, 8] * 2, 2),
        ("plain 4in16 b4", 4, [1, 2, 4, 8] * 4, 4),
        ("plain 3in15 b5", 5, [1, 2, 4, 8, 16] * 3, 3),
        ("1in4 (7,11,13,14) b4", 4, [7, 11, 13, 14], 1),
        ("a7 2in8 two reps b4", 4, [1, 7, 2, 11, 4, 13, 8, 14], 2),
        ("a8 4in16 2x2 reps b4", 4, [1, 7, 2, 11, 4, 13, 8, 14] * 2, 4),
        ("U3 2in7 b3", 3, list(range(1, 8)), 2),
        ("U3^(2) 4in14 b3", 3, list(range(1, 8)) * 2, 4),
        ("U4 4in15 b4", 4, list(range(1, 16)), 4),
        ("U4^(2) 8in30 b4", 4, list(range(1, 16)) * 2, 8),
        ("U4 + e1: 4in16 b4", 4, list(range(1, 16)) + [1], 4),
        ("U4 + e1 + e2: 4in17 b4", 4, list(range(1, 16)) + [1, 2], 4),
        ("U5 8in31 b5", 5, list(range(1, 32)), 8),
        ("2in5 b3", 3, [1, 2, 4, 3, 7], 2),
    ]
    for lab, m, S, W in cases:
        clique_test(m, S, W, label=lab)
