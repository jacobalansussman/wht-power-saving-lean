#!/usr/bin/env python3
"""cx-medium / price_t74_p1.py (own code, stdlib; python3 -I).  The hub's record words (T74 line, #1839: gauge44 and
extra; per-call profiles as read by an-hub from the words' `report`, an-hub/NOTE.md sec. 1) in the template shrunk to
one pair (block 25), every helper piece on banks with IDEAL packing (no pads; widths then do not matter).
Helper levels per call = profile rank - 224*8 (data) - 64 (copies).  Fold price at m = 25 (an-ours anatomy_ours.price)."""
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../anatomy-our-lean-words')))
from anatomy_ours import price, rank, narrow
WORDS = {'T74 gauge44': {1: 2348, 2: 500, 3: 416, 4: 130, 5: 82, 8: 8},
         'T74 extra':   {1: 2340, 2: 508, 3: 420, 4: 126, 5: 82, 8: 8},
         'U8X0 (check)': {1: 2548, 2: 338, 3: 446, 4: 162, 5: 72, 8: 8}}
M, T = 25, 25
out = {}
for name, prof in WORDS.items():
    helper = rank(prof) - 224 * 8 - 64
    banks = 3 * T * helper / M
    assert banks == int(banks), (name, banks)
    H = Counter({r: 3 * T * n for r, n in prof.items()})
    W = 224 * T + int(banks)
    D = M * W - rank(dict(H)); assert D == 32 * T
    a = float(price(dict(H), W, m=M))
    out[name] = dict(helper_levels_per_call=helper, W=W, banks=int(banks), D=D, a=a)
    print(f"{name:14s}: helper levels per call {helper}, T={T}, W={W}, banks={int(banks)}, D={D}, price {a:.7e} (x{a/1.6797550e-3:.3f} of the hub)")
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'price_t74_p1.json'), 'w'), indent=1)
