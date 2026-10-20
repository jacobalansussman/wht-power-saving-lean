#!/usr/bin/env python3
"""s5_table.py (t2-block3): which symmetric block-3 shapes W in M survive the per-character screens (PROOF.md S1-S5).

Group H of order n acting through letters (Levi lifts for fixed arrays).  Shape data: a regular array orbits, b fixed
arrays (W = n a + b); j free move orbits, s fixed moves (M = n j + s).  Per character k (0 = trivial):
  hook (S4):      80 j + sigma_{-k} s >= 224 a + d_k b
  depth (S5):     k = 0:  j + s >= a + b + 2       (fixed-array pieces of the trivial character are frame-visible)
                  k != 0: j + s >= a + 2 if a >= 1, else j + s >= 3 if b >= 1
  (S1) b = s = 0: the shape is equivalent to the smaller unit a in j (reported, not counted as new).
Splits d (c'-eigenvalue zeta^k) and sigma (c-eigenvalue zeta^k) from out_s2_split.txt / out_s4_lifts.txt;
character k pairs d_k with sigma_{-k} (Q'_chi K = K P_chi).  Only exact arithmetic on integers.
"""
import math

GROUPS = {
    "Z7 Singer (Levi)":            (7, [32] * 7, None),                 # no fixed moves at block 3 (S3)
    "Z3 3-cycle, lift A":          (3, [76, 74, 74], [28, 26, 26]),
    "Z3 3-cycle, lift B":          (3, [74, 74, 76], [26, 28, 26]),
    "Z2 transposition, lift A":    (2, [104, 120], [36, 44]),
    "Z2 transposition, lift B":    (2, [112, 112], [40, 40]),
    "Z2 transposition, lift C":    (2, [120, 104], [44, 36]),
    "Z2 Pauli X^v":                (2, [112, 112], [40, 40]),
}

def survivors(W, M, n, d, sig):
    out = []
    for a in range(W // n + 1):
        b = W - n * a
        for j in range(M // n + 1):
            s = M - n * j
            if sig is None and s:
                continue
            ok = True
            for k in range(n):
                cap = 80 * j + (sig[(-k) % n] * s if sig else 0)   # fixed registers enter character k as Q'_k: c-eigenvalue zeta^-k
                if cap < 224 * a + d[k] * b:
                    ok = False; break
            if not ok:
                continue
            if j + s < a + b + 2:
                continue
            if a >= 1 and j + s < a + 2:
                continue
            if a == 0 and b >= 1 and j + s < 3:
                continue
            out.append((a, b, j, s, "equiv. smaller unit %d in %d" % (a, j) if (b == 0 and s == 0) else "new"))
    return out

lines = []
def say(x):
    print(x); lines.append(x)

targets = [(7, 21), (2, 5), (5, 14), (6, 17), (7, 20), (10, 28), (10, 29), (12, 35), (15, 42), (15, 44), (17, 49), (21, 62), (22, 63),
           (24, 70), (25, 70), (35, 98)]
for (W, M) in targets:
    sv = ln = math.log(3 * W / M) / math.log(3)
    say("== %d in %d (saving %.4f; hook needs M >= %.1f)" % (W, M, sv, 2.8 * W))
    for gname, (n, d, sig) in GROUPS.items():
        S = survivors(W, M, n, d, sig)
        new = [x for x in S if x[4] == "new"]
        eq = [x for x in S if x[4] != "new"]
        if not S:
            say("   %-28s EXCLUDED (no (a,b,j,s) passes)" % gname)
        else:
            say("   %-28s %d new shapes%s%s" % (gname, len(new),
                (": " + ", ".join("(a%d b%d j%d s%d)" % x[:4] for x in new[:6]) + (" ..." if len(new) > 6 else "")) if new else "",
                ("; " + eq[0][4]) if eq else ""))

# all one-move-saving shapes W in 3W-1, W <= 40, that survive with Z7 (the only group whose move orbits are all free)
say("\nZ7 Singer survivors among W in M, 2.8W <= M <= 3W-1, W <= 40 (new shapes only):")
for W in range(5, 41):
    for M in range(math.ceil(2.8 * W - 1e-9), 3 * W):
        new = [x for x in survivors(W, M, 7, [32] * 7, None) if x[4] == "new"]
        if new:
            say("   %d in %d (saving %.4f): %s" % (W, M, math.log(3 * W / M) / math.log(3), ", ".join("(a%d b%d j%d)" % x[:3] for x in new)))
open("out_s5_table.txt", "w").write("\n".join(lines) + "\n")
