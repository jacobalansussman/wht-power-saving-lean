#!/usr/bin/env python3
"""t_lifts.py (rm2-phase): refusal calibration of the frame-family test in Lemma N's block, exact, from the rows only.
 block 2: search a lift (Y_U spanned by integer rows reducing mod 2 to a basis of U) of the 15 frames with excess 0;
          then abar_2/abar_1 = 2, so C[Sp(4,2)] has no unit below 2W: every saving at block 2 is refused.
 block 3: w-bounds' lift cert_E162 (135 frames of the seven-point model, 3 integer rows of length 7 each), re-checked
          here with my own code: rows mod 2 span a 3-dim totally isotropic subspace (form |A cap B| mod 2 on even
          subsets of 7 points, a nondegenerate alternating form, so these are the 135 Lagrangians of a symplectic
          F_2^6); distance from mod-2 intersections; delta from exact integer ranks.  Ratio abar_3/abar_1 > 5/2 refuses
          engine 2 in 5 in C[Sp(6,2)] (one transvection per paid move)."""
import json, random, sys
from fractions import Fraction as Fr
from spg import Sp
from t_avg_lib import rank_int

def f2_span(vs):
    S = {0}
    for v in vs:
        if v not in S: S |= {a ^ v for a in S}
    return frozenset(S)

# ---------------- block 2
sp = Sp(2); rng = random.Random(3)
from spg import lagrangians
frames, _ = lagrangians(sp)
def bits(v, n): return [(v >> i) & 1 for i in range(n)]
bases = {}
for U in frames:
    nz = sorted(u for u in U if u)
    bases[U] = [(a, b) for a in nz for b in nz if a < b]
adj = [(i, j) for i in range(len(frames)) for j in range(i + 1, len(frames)) if sp.dist(frames[i], frames[j]) == 1]
def lift(U, choice):
    (a, b), sa, sb = choice
    ra = [x * s for x, s in zip(bits(a, 4), sa)]; rb = [x * s for x, s in zip(bits(b, 4), sb)]
    return [ra, rb]
def rand_choice(U):
    return (rng.choice(bases[U]), [rng.choice([1, -1]) for _ in range(4)], [rng.choice([1, -1]) for _ in range(4)])
best = None
for restart in range(200):
    ch = [rand_choice(U) for U in frames]
    rows = [lift(U, c) for U, c in zip(frames, ch)]
    def excess(i, j): return rank_int(rows[i] + rows[j]) - 2 - 1
    E = sum(excess(i, j) for i, j in adj)
    for it in range(400):
        if E == 0: break
        k = rng.randrange(len(frames)); old = (ch[k], rows[k])
        ch[k] = rand_choice(frames[k]); rows[k] = lift(frames[k], ch[k])
        E2 = sum(excess(i, j) for i, j in adj)
        if E2 <= E: E = E2
        else: ch[k], rows[k] = old
    if best is None or E < best[0]: best = (E, [r[:] for r in rows])
    if E == 0: break
E, rows = best
# exact profile
ok_red = all(f2_span([sum((abs(x) % 2) << i for i, x in enumerate(r)) for r in R]) == U for R, U in zip(rows, frames))
prof = {}
for i in range(len(frames)):
    for j in range(len(frames)):
        if i == j: continue
        d = sp.dist(frames[i], frames[j]); prof.setdefault(d, []).append(rank_int(rows[i] + rows[j]) - 2)
ab = {d: Fr(sum(v), len(v)) for d, v in prof.items()}
print(f"block 2: lift found with total edge excess {E} (restart {restart}); rows reduce to their frames: {ok_red}; "
      f"abar_1 = {ab[1]}, abar_2 = {ab[2]}; ratio {ab[2] / ab[1]}")
if ab[2] / ab[1] == 2:
    print("   => in C[Sp(4,2)] every unit has M >= 2W: no saving at block 2 (agrees with the Lean-checked engine bound).")
with open("lift_m2_perfect.json", "w") as f:
    json.dump({"frames_pauli_labels": [sorted(U) for U in frames], "rows": rows}, f)

# ---------------- block 3: the E162 file
src = sys.argv[1] if len(sys.argv) > 1 else "../../../../notes/network-coding/cert_E162_bitmodel.json"
d = json.load(open(src)); R = d["rows"]
def red(r): return sum((abs(x) % 2) << i for i, x in enumerate(r))
def om7(a, b): return bin(a & b).count("1") & 1
spans = []; ok = True
for rows3 in R:
    vs = [red(r) for r in rows3]
    S = f2_span(vs)
    ok &= (len(S) == 8) and all(bin(v).count("1") % 2 == 0 for v in vs) and all(om7(a, b) == 0 for a in vs for b in vs)
    spans.append(S)
ok &= len(set(spans)) == 135
dimY = [rank_int(r) for r in R]
prof = {}
for i in range(135):
    for j in range(135):
        if i == j: continue
        k = len(spans[i] & spans[j]).bit_length() - 1; dd = 3 - k
        prof.setdefault(dd, {}); v = rank_int(R[i] + R[j]) - dimY[j]
        prof[dd][v] = prof[dd].get(v, 0) + 1
ab = {dd: Fr(sum(k * c for k, c in h.items()), sum(h.values())) for dd, h in prof.items()}
print(f"block 3: E162 file: 135 distinct Lagrangians of the seven-point symplectic space, rows reduce correctly: {ok}; "
      f"dims {sorted(set(dimY))}")
for dd in sorted(prof): print(f"   distance {dd}: delta histogram {dict(sorted(prof[dd].items()))}; abar = {ab[dd]}")
r = ab[3] / ab[1]
print(f"   ratio abar_3/abar_1 = {r} = {float(r):.6f}; 2 in 5 needs 5*abar_1 >= 2*abar_3: {5 * ab[1]} >= {2 * ab[3]}? "
      f"{5 * ab[1] >= 2 * ab[3]}  => engine 2 in 5 REFUSED in C[Sp(6,2)]: {5 * ab[1] < 2 * ab[3]}")
