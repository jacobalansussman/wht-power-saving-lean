"""xc3-cert3: EXACT check of a lift family on the 891 alpha-frames (block 3), over Q(omega).

Input: lift_in.txt + lift_in_base.json (prep_lift.py: Fermat lifts for the 405 Fermat frames, Teichmueller lifts of the
F4-RREF for the others).  For every frame U the lift Y_U has the F4-RREF pivot block = I and non-pivot entries in
Z[omega]_(2) congruent to the F4 entries mod 2, so Y_U reduces to U, and for every pair rank_K [Y_U; Y_V] >=
rank_F4 [U; V] = 3 + dist(U, V): delta >= dist, delta = 3 on opposite pairs.  This script computes rank_K exactly
(fractions in Q(omega)) on all 18711 edges, giving E = sum (delta - 1) and the bound 3 / (1 + E/18711) (Averaging
Lemma over Gamma, transitive on ordered pairs at each distance: geom891.py).
"""
import json, sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def qmul(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def qsub(x, y): return (x[0] - y[0], x[1] - y[1])
def qinv(x):
    a, b = x; n = a * a - a * b + b * b
    return ((a - b) / n, -b / n)
def qzero(x): return x[0] == 0 and x[1] == 0

def rank_q(M):
    M = [list(r) for r in M]; R = len(M); C = len(M[0]); r = 0
    for c in range(C):
        pv = next((i for i in range(r, R) if not qzero(M[i][c])), None)
        if pv is None: continue
        M[r], M[pv] = M[pv], M[r]
        inv = qinv(M[r][c])
        for i in range(r + 1, R):
            if not qzero(M[i][c]):
                f = qmul(M[i][c], inv)
                M[i] = [qsub(x, qmul(f, y)) for x, y in zip(M[i], M[r])]
        r += 1
    return r

if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "lift_in"
    base = json.load(open(tag + "_base.json"))
    lines = open(tag + ".txt").read().split("\n")
    P, OM, NF, NE = map(int, lines[0].split())
    edges = [tuple(map(int, l.split())) for l in lines[1 + NF:1 + NF + NE]]
    rows = {}
    for f in range(NF):
        b = base[str(f)]
        R = []
        for r in range(3):
            row = [(Fr(0), Fr(0))] * 6
            row[b["piv"][r]] = (Fr(1), Fr(0))
            for k, c in enumerate(b["npv"]):
                a1, b1 = b["base"][r][k]; row[c] = (Fr(a1), Fr(b1))
            R.append(row)
        rows[f] = R
    from collections import Counter
    H = Counter()
    for (u, v) in edges:
        H[rank_q(rows[u] + rows[v]) - 4] += 1
    E = sum(x * c for x, c in H.items())
    assert min(H) >= 0
    print(f"edges {len(edges)}; excess histogram {dict(sorted(H.items()))}; E = {E}")
    print(f"abar_1 = 1 + E/18711 = {Fr(18711 + E, 18711)} ; abar_3 = 3 ; moves >= 3W/abar_1 = {Fr(3 * 18711, 18711 + E)} W"
          f" = {3 * 18711 / (18711 + E):.6f} W  [exact over Q(omega)]")
