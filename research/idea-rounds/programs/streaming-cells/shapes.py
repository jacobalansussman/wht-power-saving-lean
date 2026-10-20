#!/usr/bin/env python3
"""Arithmetic only, from upheld hook values (RESULTS.md T-section, L7, L9, N1):
block 3 hook: cap 80; demand of a single-term target at distance 1, 2, 3 = 80, 144, 224.
block 4 hook (k = 8): cap 27,456; demand at distance 1, 2, 3, 4 = 27,456, 48,576, 68,320 (computed), >= 86,592.
A minimal streaming cell with delays in {0,1} is a block program whose arrays are owed single-term stages
(proved in NOTE.md section 3); Theorem R then needs c * cap >= sum of stage demands.
Also: unknowns versus equations for cells in the no-phase game (question 1)."""
from math import log
from itertools import combinations_with_replacement as cwr

DEM = {3: {1: 80, 2: 144, 3: 224}, 4: {1: 27456, 2: 48576, 3: 68320, 4: 86592}}
CAP = {3: 80, 4: 27456}

def chainable(ranks, m):
    """Can the stage lengths be split into groups each summing to m (chains composing to K)?"""
    rs = sorted(ranks, reverse=True)
    if sum(rs) % m:
        return False
    def place(i, bins):
        if i == len(rs):
            return all(b == m for b in bins)
        for k in range(len(bins)):
            if bins[k] + rs[i] <= m:
                bins[k] += rs[i]
                if place(i + 1, bins):
                    return True
                bins[k] -= rs[i]
            if bins[k] == 0:
                break
        return False
    return place(0, [0] * (sum(rs) // m))

def tau(ranks, M):
    lo, hi = 0.0, 1.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if sum(k ** mid for k in ranks) > M:
            hi = mid
        else:
            lo = mid
    return lo

print("Mixed / split-stage shapes one move below plain that the hook does NOT exclude")
print("(ranks = per-array stage lengths; 'stream' = saving if used as a streaming cell with full")
print(" transforms assembled from the stages; 'mixed' = 1 - tau of the tau-theorem M18)")
for m in (3, 4):
    found = []
    for W in range(1, 5):
        for ranks in cwr(range(1, m + 1), W):
            plain = sum(ranks)
            c = plain - 1
            need = sum(DEM[m][k] for k in ranks) / CAP[m]
            if c >= need and any(k >= 2 for k in ranks):
                slack = c - need
                stream = log(plain / c) / log(m) if chainable(ranks, m) else float('nan')
                found.append((W, ranks, c, slack, 1 - tau(ranks, c), stream))
    print(f"  block {m}: {len(found)} shapes with at most 4 arrays")
    for W, ranks, c, slack, mixed, stream in sorted(found, key=lambda t: (t[2], -t[4]))[:14]:
        s = f"{stream:.4f}" if stream == stream else "  -   "
        print(f"    ranks {ranks} in {c} moves: hook slack {slack:.3f} moves; mixed value {mixed:.4f}; stream {s}")

print()
print("Question 1: unknowns versus equations, no-phase game, cyclic address group of order n,")
print("cell with p inputs, r states, c moves, generic C[H] coefficients; equations = forms x Markov lags")
rows = [(4, 15, 2, 0, 7, 0), (4, 15, 2, 1, 7, 1), (4, 15, 2, 2, 7, 1), (5, 31, 1, 1, 4, 1), (5, 31, 1, 2, 4, 2),
        (6, 63, 1, 2, 5, 2), (7, 127, 1, 2, 6, 2), (8, 255, 1, 3, 7, 3)]
for m, n, p, r, c, d in rows:
    U = n * (c * (p + r) + c * (c - 1) // 2 + (p + r) * (p + r + c))
    forms = 2 ** (m * (m + 1) // 2)
    E = forms * n * p * p * (d + n * r + 1)
    print(f"  block {m}, n={n}, p={p}, r={r}, c={c}, delay={d}: unknowns {U:,}  scalar equations ~{E:.2e}"
          f"  (forms {forms:.2e})")
