# e3-boost lp3.py: boosters with small ranks on LARGE registers (spare rows used only by free steps): registers
# m = 6..12, ranks <= rmax.  Same LP as lp2.boost_min; certificates: rigorous + argued (projections).
import sys
import certlib
from lp import rows
from lp2 import boost_min
grid = [i / 100 for i in range(1, 100)] + [0.9983]
st = certlib.RIGOROUS + ("argued",)
for m in range(6, 13):
    R, names = rows(m, st)
    for rmax in (4, 5, 6):
        if rmax > m: continue
        neg = [t for t in grid if boost_min(m, t, R, rmax)[0] < -1e-9]
        print(f"register {m}, ranks <= {rmax}: {len(names)} certificate rows; boosters negative at {len(neg)}/{len(grid)} "
              f"grid points" + (f", from {min(neg)} to {max(neg)}" if neg else ""), flush=True)
