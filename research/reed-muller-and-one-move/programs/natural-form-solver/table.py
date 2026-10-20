"""table.py: summary of all 8/9/10-in-31 random-start runs (run.py logs): per target set, runs, best residual on
all forms, most arrays served, plain packing number."""
import sys, re, glob, itertools
from packing2 import packing
fam = {}
for fn in sorted(glob.glob("jobs/*/remote-out/*/remote.log")):
    txt = open(fn).read().replace("\n ", " ")
    ks = None; slice_ = False
    for line in txt.splitlines():
        m = re.search(r"\] m5 ks=\[([0-9, ]+)\] dirs=all", line)
        if m: ks = tuple(int(v) for v in m.group(1).split(",")); slice_ = False; continue
        m = re.search(r"\] m5 slice I x(\d+) dirs=all", line)
        if m: ks = tuple([0] * int(m.group(1))); slice_ = True; continue
        m = re.search(r"\] m5 .*dirs=frames", line)
        if m: ks = None; continue
        m = re.search(r"random survivor ks=\[([0-9, ]+)\]", line)
        if m: ks = tuple(int(v) for v in m.group(1).split(",")); continue
        m = re.search(r"\] RESTART \d+ subset \S+ ALL (\S+) per-array \[([^\]]+)\]", line)
        if m and ks is not None:
            pa = [float(v) for v in m.group(2).split()]
            if len(pa) != len(ks): continue
            key = (ks, slice_)
            r = fam.setdefault(key, [0, 1e99, 0])
            r[0] += 1; r[1] = min(r[1], float(m.group(1))); r[2] = max(r[2], sum(v < 1e-10 for v in pa))
tot = 0
for (ks, sl), (n, best, srv) in sorted(fam.items(), key=lambda kv: (len(kv[0][0]), kv[1][1])):
    pn = packing(list(ks))[0]
    tot += n
    print(f"| {len(ks)} in 31 | {'slice, all I' if sl else ','.join(map(str, ks))} | {n} | {best:.3e} | {srv} | {pn} |")
print("total runs", tot)
