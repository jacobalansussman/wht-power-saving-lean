# v-e3-boost: e3-boost's lp3 (small ranks on registers 6..12) re-run with RIGOROUS certificates only (no argued
# projections), to see whether the NOTE's "[computed]" label for "all jobs <= 5 rows, registers <= 12" needs the
# projections.  Imports e3-boost's certlib/lp/lp2 read-only.
import sys
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
import certlib
from lp import rows
from lp2 import boost_min
grid = [i / 50 for i in range(1, 50)] + [0.9983]
for m in range(6, 13):
    R, names = rows(m, certlib.RIGOROUS)
    for rmax in (4, 5):
        neg = [t for t in grid if boost_min(m, t, R, rmax)[0] < -1e-9]
        print(f"RIGOROUS only: register {m}, ranks <= {rmax}: {len(names)} rows; booster-negative at {len(neg)}/{len(grid)}"
              + (f" ({min(neg)}..{max(neg)})" if neg else ""), flush=True)
