"""summarize.py -- parse the RESULT lines of the Python-call timing files into the report table.
For each cell: the SME program (fc_sme), the best baseline among numpy / scipy (auto, oa, fft) / vDSP, and the change.
Usage: python summarize.py FILE... ; prints markdown rows (best of 9; medians in brackets)."""
import re
import sys
from collections import defaultdict

BASE = ("np", "sp_auto", "sp_oa", "sp_fft", "vdsp")
cells = defaultdict(dict)
order = []
for fn in sys.argv[1:]:
    for line in open(fn):
        if not line.startswith("RESULT"):
            continue
        m = re.match(r"RESULT N=(\d+) M=(\d+) (\w+) (\w+) +(\w+) +best +([\d.]+) median +([\d.]+)", line)
        if not m:
            continue
        N, M, mode, dt, prog, best, med = m.groups()
        key = (dt, mode, int(N), int(M))
        if key not in cells:
            order.append(key)
        cells[key][prog] = (float(best), float(med), fn)
order.sort()
print("| dtype | mode | N | M | fc64 SME | best baseline | baseline ns | change best (median) | numpy | other fc64 |")
print("|---|---|---|---|---|---|---|---|---|---|")
for key in order:
    c = cells[key]
    if "fc_sme" not in c and "fc_auto" not in c:
        continue
    me = "fc_sme" if "fc_sme" in c else "fc_auto"
    bl = [(c[b][0], b) for b in BASE if b in c]
    if not bl:
        continue
    bb, bn = min(bl)
    bmed = c[bn][1]
    s, smed = c[me][0], c[me][1]
    ch = 100 * (s / bb - 1)
    chm = 100 * (smed / bmed - 1)
    other = ", ".join(f"{p} {c[p][0]:.3f}" for p in c if p.startswith("fc_") and p != me)
    npv = f"{c['np'][0]:.3f}" if "np" in c else "-"
    print(f"| {key[0]} | {key[1]} | {key[2]:,} | {key[3]} | {s:.3f} ({smed:.3f}) | {bn} | {bb:.3f} ({bmed:.3f}) | "
          f"{ch:+.0f} % ({chm:+.0f} %) | {npv} | {other} |")
