#!/usr/bin/env python3
"""t_scalar.py (p6-prove): ledger.py's normalisation modulo {1,i,-1,-i} is not the every-f model.

Word (W = 1, block 1, registers c, s, t): s := c; move s (s = D1 x); c := c + s (c = x + D1 x);
t := c0-copy of x made at the start, then XZXZ on t (as matrices XZXZ = -I, a GROUP element acting as (-1)^f);
c := c + t.  ledger.py folds -I into the coefficient and reports the unit D1 x; in the true model the content is
x + (-I)x + D1 x, three distinct group elements, not a unit (at even f the first two add up to 2x).
"""
import sys
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "ledger"))
from ledger import Engine, run
from checks import true_unit
E = Engine(1)
Z = None
# Z = (P1_0)^2 up to a scalar; X1 = X
w = [("mix", 2, [(0, 1, None)]), ("mix", 1, [(0, 1, None)]), ("move", 1, 1), ("mix", 0, [(0, 1, None), (1, 1, None)])]
w += [("free", 2, "X1"), ("free", 2, "P1_0"), ("free", 2, "P1_0"), ("free", 2, "X1"), ("free", 2, "P1_0"), ("free", 2, "P1_0")]
w += [("mix", 0, [(0, 1, None), (2, 1, None)])]
res = run(E, 1, 3, w)
print("ledger.py unit verdict:", res["unit"], " contents of c:", res["contents"][0])
print("true every-f model unit verdict:", true_unit(E, 1, 3, w))
