#!/usr/bin/env python3
"""s12_s3.py (t2-block3, stdlib): S3 = the coordinate permutations (an exact symmetry of K and of the letters) acting on
W in M at block 3, per irreducible representation of S3 (PROOF.md S10).

Array orbit with stabiliser S contributes dem_rho = <rho|_S, (N/E_N)|_S>_S, move orbit with stabiliser T (letter m fixed
by T) contributes cap_rho = <rho|_T, ((N + mN)/N)|_T>_T (Frobenius reciprocity), from the splits of s2/s4:
  N/E_N:      3-cycle (76,74,74) -> trace 2;  transposition (104,120) -> trace -16
  (N+mN)/N:   3-cycle (28,26,26) -> trace 2;  transposition (36,44)  -> trace -8   (dir7; all tested fixed letters agree)
Screens: per irrep rho in {1, sign, 2-dim}: sum caps >= sum demands; orbit count (S5') for S3 and for its subgroups
Z3 and Z2: (#move orbits) >= (#array orbits) + 2.
"""
import math, itertools
# per orbit type: (size, dem(1, sgn, std))
ARR = {"S3": (1, (30, 46, 74)), "Z3": (2, (76, 76, 148)), "Z2": (3, (104, 120, 224)), "1": (6, (224, 224, 448))}
MOV = {"S3": (1, (10, 18, 26)), "Z3": (2, (28, 28, 52)), "Z2": (3, (36, 44, 80)), "1": (6, (80, 80, 160))}
# sanity: dims
for k, (sz, d) in ARR.items():
    assert d[0] + d[1] + 2 * d[2] == 224 * sz
for k, (sz, d) in MOV.items():
    assert d[0] + d[1] + 2 * d[2] == 80 * sz
# orbits of subgroups inside an S3-orbit with stabiliser S: (#Z3-orbits, #Z2-orbits) for the subgroup <3-cycle> and <one transposition>
SUB = {"S3": (1, 1), "Z3": (2, 1), "Z2": (1, 2), "1": (2, 3)}
types = ["S3", "Z3", "Z2", "1"]

def decomps(total, table):
    sizes = [table[t][0] for t in types]
    res = []
    def rec(i, left, cur):
        if i == len(types):
            if left == 0: res.append(tuple(cur))
            return
        for c in range(left // sizes[i] + 1):
            rec(i + 1, left - c * sizes[i], cur + [c])
    rec(0, total, [])
    return res

lines = []
def say(x):
    print(x); lines.append(x)
for (W, M) in [(6, 18), (2, 5), (5, 14), (6, 17), (7, 20), (10, 28), (10, 29), (12, 35), (15, 42), (15, 44), (21, 62)]:
    found = []
    for A in decomps(W, ARR):
        for Mv in decomps(M, MOV):
            ok = True
            for r in range(3):
                cap = sum(Mv[i] * MOV[t][1][r] for i, t in enumerate(types))
                dem = sum(A[i] * ARR[t][1][r] for i, t in enumerate(types))
                if cap < dem: ok = False; break
            if not ok: continue
            if sum(Mv) < sum(A) + 2: continue
            for s in range(2):   # subgroup Z3 (s=0) and Z2 (s=1)
                if sum(Mv[i] * SUB[t][s] for i, t in enumerate(types)) < sum(A[i] * SUB[t][s] for i, t in enumerate(types)) + 2:
                    ok = False
            if ok:
                found.append((A, Mv))
    tag = lambda v: "+".join("%dx%s" % (c, {"S3": "1", "Z3": "2", "Z2": "3", "1": "6"}[t]) for c, t in zip(v, types) if c)
    say("%d in %d (saving %.4f): %s" % (W, M, math.log(3 * W / M) / math.log(3),
        ("%d S3-shapes; e.g. " % len(found) + "; ".join("arrays %s | moves %s" % (tag(a), tag(m)) for a, m in found[:4])) if found else "EXCLUDED for S3"))
open("out_s12_s3.txt", "w").write("\n".join(lines) + "\n")
