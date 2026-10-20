# e3-boost projall.py: carry EVERY symmetric family (A, B, C_k, distance sets) and every projected profile from
# blocks 3..7 up to block 8 (status 'argued2'), then rerun the booster LP on registers 6..8 (ranks <= rmax).
import sys
import certlib as C
from lp import rows
from lp2 import boost_min, cover
base = [(n, M, P, s) for (n, M, P, s) in list(C.PROFILES) if 3 <= M <= 7 and s != "sampled" and s != "conj"
        and not n.startswith("P")]
cnt = 0
for (n, M, P, s) in base:
    cur = list(P)
    for M2 in range(M + 1, 9):
        cur = C.project(cur)
        if max(cur) > 0: C.add(f"Q{n}to{M2}", M2, cur, "argued2"); cnt += 1
print("projected profiles added:", cnt)
st = C.RIGOROUS + ("argued", "argued2")
grid = [i / 50 for i in range(1, 50)] + [0.9983]
for m in (6, 7, 8):
    R, names = rows(m, st)
    for rmax in range(5, m + 1):
        neg = [t for t in grid if boost_min(m, t, R, rmax)[0] < -1e-9]
        print(f"register {m}, ranks <= {rmax}: {len(names)} rows; boosters negative at {len(neg)}/{len(grid)}"
              + (f", from {min(neg)} to {max(neg)}" if neg else ""), flush=True)
    e, v = cover(m, R, 0.9)
    print(f"   cover [0.9, 0.9983]: eps = {e:.3g}, v = {dict((d + 1, round(float(x), 5)) for d, x in enumerate(v) if abs(x) > 1e-7)}")
