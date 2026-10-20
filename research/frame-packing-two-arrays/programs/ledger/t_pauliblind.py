#!/usr/bin/env python3
"""t_pauliblind.py: a sub-class (S2) example -- frame-coherent but strictly incoherent (Pauli-blind packing):
r0 := (1 + Z^v) y0, one move, then the twin is removed with a scratch copy.  Wasteful by one move.  Shows the frame
ledger sees no packing while the strict ledger does; both balance."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger import *
for m in (2, 3, 4):
    E = Engine(m); v = (1 << m) - 1
    w = [("mix", 1, [(0, 1, None)]), ("free", 1, "P%d_0" % v), ("free", 1, "P%d_0" % v),      # s = Z^v y0 (up to sign)
         ("mix", 0, [(0, 1, None), (1, 1, None)]),                                               # r0 = (1 + Z^v) y0
         ("move", 0, 1), ("move", 1, 1),                                                         # packed move; twin's own move
         ("mix", 0, [(0, 1, None), (1, -1, None)])]                                              # r0 = D1 y0
    w += [("move", 0, 1 << j) for j in range(1, m)]
    st = run(E, 1, 2, w, start={0: 0}, deliver={0: 0}, strict=True)
    fr = run(E, 1, 2, w, start={0: 0}, deliver={0: 0})
    gs, why_s = find_strict_gauge(E, 1, 2, w, start={0: 0})
    gf, why_f = find_coherent_gauge(E, 1, 2, w, start={0: 0})
    print("m=%d M=%d unit=%s | frame: %s, P=%d B=%d ok=%s | strict: %s, P=%d B=%d ok=%s"
          % (m, st["M"], st["unit"], why_f, fr["P"], fr["B"], fr["unit_ok"], why_s, st["P"], st["B"], st["unit_ok"]))
