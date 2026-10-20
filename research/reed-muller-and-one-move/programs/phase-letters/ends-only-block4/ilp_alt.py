"""ilp_alt.py (rm2-phase, restart session): rm2-exact's ILP piece + term-rank bound (rc_ilp.py, copied unchanged) on
classes with targets of any kind; the verifier's second parity is replaced by the decomposition-based one of
tile_alt.py (alternating targets need m + 1 rank-ones), everything else is rc_ilp's own code.
usage: ilp_alt.py m pass_lines.txt [which=all|even] [maxforms]
"""
import sys, re, time
import numpy as np
import rc_ilp
from xbase import Field, Forms, System
from tile_alt import target_matrices, second_parity as sp2
rc_ilp.second_parity = lambda Fd, fo, q, T: sp2(fo, q, T)

m = int(sys.argv[1]); path = sys.argv[2]; which = sys.argv[3] if len(sys.argv) > 3 else "all"
maxforms = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
mats, nonalt = target_matrices(m); Fd = Field(m); fo = Forms(m)
pat = re.compile(r"class (\d+): dirs \(([\d, ]+)\) targets \(([\d, ]+?),?\) ")
for L in open(path):
    mm = pat.search(L)
    if not mm: continue
    ci = int(mm.group(1)); zs = [int(x) for x in mm.group(2).split(",")]; tt = [int(x) for x in mm.group(3).split(",")]
    counts = {}
    for z in zs: counts[z] = counts.get(z, 0) + 1
    targets = [mats[t] for t in tt]; sy = System(fo, targets, counts)
    present = sum(1 << i for i, z in enumerate(fo.dirs) if counts.get(z, 0) > 0)
    if which == "even": cand = [F for F in range(1, fo.nf) if sy.r[F] == 0 and fo.mask[F] & present]
    else: cand = [F for F in range(1, fo.nf) if sy.s[F] > 0 and fo.mask[F] & present]
    cand.sort(key=lambda F: (-(sy.s[F] / max(1, sy.d[F])), sy.d[F])); cand = cand[:maxforms]
    t0 = time.time(); nviol = 0; hist = {}; first = None
    for F in cand:
        val, npc, info = rc_ilp.solve(fo, sy, F, 10 ** 9)
        if val is None: hist["fail"] = hist.get("fail", 0) + 1; continue
        mg = val - int(sy.s[F]); hist[mg] = hist.get(mg, 0) + 1
        if val < sy.s[F]:
            chosen, rows, cols, st = info
            tot, sF = rc_ilp.verify(Fd, fo, sy, targets, F, chosen, rows, cols); nviol += 1
            if first is None: first = (F, tot, sF, len(chosen), len(rows) + len(cols))
            break
    print(f"class {ci} targets {tuple(tt)} ({sum(t >= nonalt for t in tt)} alternating): ILP on {len(cand)} forms ({which}) "
          f"in {time.time()-t0:.0f}s; VERIFIED violations {nviol}{' first F=%d bound %d < s_F %d (%d pieces + %d lines)' % first if first else ''}; "
          f"margins {sorted(hist.items(), key=str)[:6]}", flush=True)
