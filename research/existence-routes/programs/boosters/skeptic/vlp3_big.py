# v-e3-boost: e3-boost's lp3 (ranks <= 4, 5; rigorous + argued, exactly as lp3) extended to registers 13..20, which the
# NOTE does not cover but whose conclusion "any booster needs a job of at least 6 rows" would need.
import sys
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
import certlib
from lp import rows
from lp2 import boost_min
grid = [i / 50 for i in range(1, 50)] + [0.9983]
st = certlib.RIGOROUS + ("argued",)
for m in list(map(int, sys.argv[1:])) or range(13, 21):
    R, names = rows(m, st)
    for rmax in (4, 5):
        neg = [t for t in grid if boost_min(m, t, R, rmax)[0] < -1e-9]
        print(f"rigorous+argued: register {m}, ranks <= {rmax}: {len(names)} rows; booster-negative at {len(neg)}/{len(grid)}"
              + (f" ({min(neg)}..{max(neg)})" if neg else ""), flush=True)
