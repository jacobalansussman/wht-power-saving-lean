#!/usr/bin/env python3
"""derived.py RAWDIR OUTFILE : "multiply then relinearize" = EvaluateMulCt + EvaluateRelinInplace, summed per run
(medians of repetitions), then compared OFF vs ON per pair exactly like analyse.py. SEAL's benchmark has no single
multiply-and-relinearize case, so this is a derived figure (the sum of two measured cases)."""
import os, re, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyse import load_run, key_of

rawdir, outfile = sys.argv[1], sys.argv[2]
runs = {}
for fn in sorted(os.listdir(rawdir)):
    m = re.match(r"(r\d+)-(small|big)-(off|on)\.json$", fn)
    if m:
        runs.setdefault((m.group(1), m.group(2)), {})[m.group(3)] = load_run(os.path.join(rawdir, fn))[0]
rows = {}
for (tag, g), pr in sorted(runs.items()):
    if "off" not in pr or "on" not in pr:
        continue
    for scheme in ("BFV", "CKKS"):
        sums = {}
        for mode in ("off", "on"):
            for name, v in pr[mode].items():
                n, op = key_of(name)
                if op in (f"{scheme} / EvaluateMulCt", f"{scheme} / EvaluateRelinInplace"):
                    sums.setdefault((mode, n), 0.0)
                    sums[(mode, n)] += statistics.median(v)
        for (mode, n) in list(sums):
            if mode == "off" and ("on", n) in sums:
                rows.setdefault((scheme, n), []).append((sums[("off", n)], sums[("on", n)]))
out = ["| n | scheme | multiply + relinearize OFF (us) | ON (us) | change | spread across pairs | verdict |",
       "|---|---|---|---|---|---|---|"]
for (scheme, n) in sorted(rows, key=lambda k: (k[1], k[0])):
    v = rows[(scheme, n)]
    a = statistics.median(x[0] for x in v); b = statistics.median(x[1] for x in v)
    per = [x[1] / x[0] - 1 for x in v]
    lo, hi = min(per), max(per)
    verdict = "no change" if lo <= 0 <= hi else ("faster" if hi < 0 else "slower")
    out.append(f"| {n} | {scheme} | {a:,.0f} | {b:,.0f} | {(b/a-1)*100:+.1f} % | {lo*100:+.1f} to {hi*100:+.1f} % | {verdict} |")
open(outfile, "w").write("\n".join(out) + "\n")
print("\n".join(out))
