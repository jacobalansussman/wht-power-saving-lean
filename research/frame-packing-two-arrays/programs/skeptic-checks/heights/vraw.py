#!/usr/bin/env python3
"""vraw.py (v-p6-prove): raw complex arrays only, for the units of vunits.py at given (m, f) pairs.
Usage: python3 vraw.py m:f1,f2,... [m:f...]"""
import sys
import numpy as np
from veng import Eng, raw_check
from vunits import words
rng = np.random.default_rng(11)
for arg in sys.argv[1:]:
    m, fs = arg.split(":"); m = int(m); fs = [int(x) for x in fs.split(",")]
    E = Eng(m)
    for (name, W, R, w) in words(E):
        for f in fs:
            print("m=%d %-7s f=%d  max rel err %.2e" % (m, name, f, raw_check(E, W, R, w, f, rng)), flush=True)
