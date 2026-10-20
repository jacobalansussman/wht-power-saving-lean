# Exact hook gains at block 3

**What it shows.** `exact_hook.py` (pure Python, exact over Q(i)): the hook values f(d1) = 80, f(d2 d1) = 144, f(d4 d2
d1) = 224, and the second move after a free letter: 144 after a plain phase, 152 after a controlled-S, and 160 = 2 ×
80 after a Toffoli, so a free Toffoli between two moves makes the second move perfect for the hook. `t4exact.py`: a
single-fibre partial move at block 3 has exact rank 20 in the hook, exactly its share of a whole move's cap, so
partial paid moves open nothing at block 3. `rankp.c` and `wt.py` are helpers of the track's batch runs.

**Commands.** `python3 -I exact_hook.py`, `python3 -I t4exact.py` (pure Python; minutes). Outputs in
`../../data/tightness-at-block-3/`. One program each, not repeated by anyone.

**Re-run for this revision:** no.
