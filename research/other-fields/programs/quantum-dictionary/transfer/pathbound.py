"""x-transfer, 13 Oct 2026 (python3 -I).  Path bound in the engine's own matrices, m = 1, 2:
for every g in G_m, (moves needed before) + (moves needed after) = rankB(g) + rankB(K g^-1) = m + e(L(g)),
where e = excess of the frame of g.  Expected histogram of e: |B_m| times the frame counts by excess
(m=1: e0 2 frames, e1 1; m=2: e0 5, e1 6, e2 4; see frames_geom.out section B)."""
import sys; sys.path.insert(0, __import__('os').path.dirname(__file__))
from segnorm import free_gens, on_qubit, SX, closure, eye, mul, dagger, rankB
for m in (1, 2):
    n = 1 << m
    F = free_gens(m); dirs = [on_qubit(m, j, SX) for j in range(m)]
    G = closure([eye(n)], F + dirs, m)
    K = eye(n)
    for d in dirs: K = mul(d, K)
    hist = {}
    for g in G:
        e = rankB(g, m) + rankB(mul(K, dagger(g)), m) - m
        hist[e] = hist.get(e, 0) + 1
    print(f"m={m}: |G|={len(G)}; histogram of rankB(g)+rankB(K g^-1)-m: {dict(sorted(hist.items()))}; min {min(hist)} (never below 0)")
