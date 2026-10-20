#!/usr/bin/env python3
"""p-wgen (wht45/lever, 19 Oct 2026): the residual list and local histogram of a gcert/1w word, in the `.info.json`
form k-word's pack_ceil.py reads (pad-minimal dedicated bank table, no source copies).  Residual of slot q (slot
order) exactly as k-combine's make_sched_x.py / Lean's res_ok: dim(final) - dim(start) for an exit slot, h - dim(start)
otherwise; checked against the gcert/1w `residuals` key.  Histogram = merged `blocks`, checked against
`hub_local_profile`.   usage: mk_info.py <word.gcert1w.json.gz> <out.info.json>"""
import sys, json, gzip
from collections import Counter
c = json.loads(gzip.open(sys.argv[1], "rt").read())
v, R, h = c["v"], c["R"], c["h"]
dim = [len(f) for f in c["frames"]]
E = set(c["exits"])
res = [(dim[c["final"][2 * v + q]] if 2 * v + q in E else h) - dim[c["start"][2 * v + q]] for q in range(R)]
assert res == c["residuals"], "residuals differ from the gcert/1w key"
hist = Counter()
for w in c["blocks"].values():
    for r, n in w.items():
        hist[int(r)] += n
assert {str(k): n for k, n in hist.items()} == c["hub_local_profile"], (hist, c["hub_local_profile"])
json.dump(dict(residuals=res, histogram={str(k): n for k, n in sorted(hist.items())}), open(sys.argv[2], "w"))
print("R", R, "residuals", sorted(Counter(res).items()), "hist", dict(sorted(hist.items())))
