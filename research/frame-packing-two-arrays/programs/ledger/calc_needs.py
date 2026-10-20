#!/usr/bin/env python3
"""calc_needs.py: what a unit with saving s needs, in the two readings of the ledger, and what the proved lower
bounds allow at each block (2W at block 2, 14W/5 at block 3, 41W/13 at block 4, 2.9447W for every m >= 4)."""
import math
floors = {2: 2.0, 3: 14 / 5, 4: 41 / 13}
def floor(m): return floors.get(m, 2.9447 if m >= 4 else m)
print("m  | proved floor M/W | max saving allowed | s=0.1: M/W, Delta/(Wm), avg crossing distance, extra labels/move "
      "| s=0.01: same")
for m in (2, 3, 4, 5, 6, 8, 12, 16, 25, 41, 73):
    f = floor(m); smax = math.log(m / f) / math.log(m) if f < m else 0.0
    row = []
    for s in (0.1, 0.01):
        mw = m ** (1 - s); frac = 1 - m ** (-s)
        row.append("M/W<=%.3f %s D/(Wm)>=%.4f cross>=%.3f pack>=%.4f" % (mw, "ALLOWED " if mw >= f else "EXCLUDED", frac, m * frac, m ** s - 1))
    print("%2d | %.4f | %.4f | %s | %s" % (m, f, smax, row[0], row[1]))
