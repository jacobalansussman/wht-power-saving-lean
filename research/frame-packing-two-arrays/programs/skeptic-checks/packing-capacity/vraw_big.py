#!/usr/bin/env python3
"""vraw_big.py (v-fp-capacity): the rebuilt named units (plain, pair trick, crossing pair, order-3 control) on raw
complex arrays at the large f the Mac cannot hold: m = 3 at f = 7, 8; m = 4 at f = 5, 6. Run on a Colab."""
import sys, json, random, time
import numpy as np
from vengine import Eng, raw_arrays
from vunits import words
out = []
for m, fs in ((3, (7, 8)), (4, (5, 6))):
    E = Eng(m)
    for (name, W, R, w, dl) in words(E, random.Random(0), 0):
        for f in fs:
            t0 = time.time()
            err = raw_arrays(E, W, R, w, dl, f, np.random.default_rng(100 + f))
            rec = dict(m=m, name=name, f=f, max_err=err, secs=round(time.time() - t0, 1))
            print(json.dumps(rec), flush=True); out.append(rec)
print("SUMMARY max_err=%.3e over %d runs" % (max(r["max_err"] for r in out), len(out)))
