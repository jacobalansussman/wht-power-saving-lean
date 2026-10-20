"""parse_runs.py LOG...: for each RESTART line of run.py logs, which arrays are served exactly (per-array < 1e-10),
whether their frames are pairwise disjoint (a plain packing), and the residual."""
import sys, re, itertools
from packing import frames, packing
rows = []
for fn in sys.argv[1:]:
    ks = None; tag = None; txt = open(fn).read().replace("\n ", " ")
    for line in txt.splitlines():
        m = re.search(r"\[(\S+)\] m5 ks=\[([0-9, ]+)\]", line)
        if m: ks = [int(v) for v in m.group(2).split(",")]; tag = m.group(1); continue
        m = re.search(r"random survivor ks=\[([0-9, ]+)\]", line)
        if m: ks = [int(v) for v in m.group(1).split(",")]; continue
        m = re.search(r"\[(\S+)\] m5 slice I x(\d+)", line)
        if m: ks = [0] * int(m.group(2)); tag = m.group(1); continue
        m = re.search(r"\[(\S+)\] RESTART (\d+) subset (\S+) ALL (\S+) per-array \[([^\]]+)\]", line)
        if m and ks is not None:
            pa = [float(v) for v in m.group(5).split()]
            served = [i for i, v in enumerate(pa) if v < 1e-10]
            fr = frames(ks)
            disj = all(not (fr[i] & fr[j]) for i, j in itertools.combinations(served, 2))
            pn = packing(ks)[0]
            rows.append((m.group(1), len(ks), float(m.group(4)), len(served), disj, pn))
            print(f"{m.group(1):14s} W={len(ks):2d} ALL {float(m.group(4)):.4e} served {served} disjoint-frames {disj} packing-number {pn} "
                  f"unserved-min {min([v for v in pa if v >= 1e-10], default=0):.1f}")
