#!/usr/bin/env python3
"""vsharp.py (v-fp-capacity): a block-2 unit (W = 1, two scratch registers, Delta = -1) on which the image-frame
count is NOT sharper than the Sp-label count: P_img = 0 > P_asc(Sp) = -1, because the class-graph liveness of
Definition 1.1 makes a coarse class live whose only fine member is cancelled. Both bounds hold (Delta = -1)."""
import numpy as np
from vengine import Eng, run, check_unit, raw_arrays, analyse, moves_report
E = Eng(2); D = E.dirm; s = E.addr([1, 0])
w = [("mix", 1, [(0, 1, None)]),                       # r1 := y
     ("move", 1, D(2)),                                # r1 = D2 y
     ("mix", 2, [(0, 1, None), (0, 1, s)]),            # r2 := y + s y
     ("move", 2, D(2)),                                # r2 = D2 y + D2 s y
     ("mix", 3, [(2, 1, None), (1, -1, None)]),        # r3 := r2 - r1 = D2 s y
     ("move", 3, D(1)),                                # r3 = D1 D2 s y = K s y
     ("free", 3, s)]                                   # r3 = s K s y = K y
W, R, dl = 1, 4, {3: 0}
snaps = run(E, W, R, w)
print("formal identity:", check_unit(E, W, snaps, dl), " raw arrays f=1..8 max err:",
      max(raw_arrays(E, W, R, w, dl, f, np.random.default_rng(f)) for f in range(1, 9)))
L0 = E.L0(); Lf = E.frame(tuple(1 << j for j in range(4)), L0)
lev = lambda F: 2 - (len(F & Lf).bit_length() - 1)
for name, cls, lf in (("Sp labels", lambda g, l: E.sp(g), lambda S: lev(E.frame(S, L0))),
                      ("image frames", lambda g, l: E.frame(E.sp(g), L0), lev)):
    nodes, live = analyse(E, W, R, w, snaps, dl, cls)
    rep = moves_report(E, W, R, w, snaps, nodes, live, cls, lf)
    per = [sum(1 for a, b in st if b == a + 1) - 1 for k, st in rep]
    print("%-13s live classes per move %s  ascents-1 per move %s  P_asc = %d   (M = 3, Delta = %d)"
          % (name, [k for k, _ in rep], per, sum(per), W * 2 - 3))
