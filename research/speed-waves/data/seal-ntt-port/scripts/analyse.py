#!/usr/bin/env python3
"""analyse.py RAWDIR OUTFILE [TAGPREFIX [A B]]

A and B are the two run modes compared (default off and on: "OFF" is A, "ON" is B; the control experiment also has
onscalar).

Reads the Google Benchmark JSON files written by run_pair.sh (TAG-GROUP-off.json / TAG-GROUP-on.json) and prints
one table: per benchmark, the OFF median, the ON median, the change, and the spread of the change across pairs.

- For each run (one file) and benchmark: the median of its repetitions' real_time (microseconds per operation;
  SEAL's benchmark runs 10 iterations per repetition).
- A pair is the OFF run and the ON run of the same TAG and GROUP, taken back to back under one timing lock.
- OFF median / ON median: the median, across pairs, of the per-run medians.
- change: ON median / OFF median - 1.
- spread: the lowest and highest per-pair change (ON/OFF - 1 within each pair).
- verdict: "no change" when the spread includes 0 (the pairs disagree on the direction); otherwise faster/slower.
Only the Python standard library is used.
"""
import json
import os
import re
import statistics
import sys

ORDER_OP = [
    "UTIL / NTTForwardLowLevelLazy", "UTIL / NTTInverseLowLevelLazy", "UTIL / NTTForwardLowLevel",
    "UTIL / NTTInverseLowLevel", "UTIL / NTTForward", "UTIL / NTTInverse",
    "BFV / EncryptSecret", "BFV / EncryptPublic", "BFV / Decrypt", "BFV / EncodeBatch", "BFV / DecodeBatch",
    "BFV / EvaluateMulCt", "BFV / EvaluateSquare", "BFV / EvaluateRelinInplace", "BFV / EvaluateRotateRows",
    "BFV / EvaluateRotateCols",
    "CKKS / EncryptSecret", "CKKS / EncryptPublic", "CKKS / Decrypt", "CKKS / EncodeDouble", "CKKS / DecodeDouble",
    "CKKS / EvaluateMulCt", "CKKS / EvaluateSquare", "CKKS / EvaluateRelinInplace", "CKKS / EvaluateRescaleInplace",
    "CKKS / EvaluateRotate", "CKKS / EvaluateAddCt",
]


def load_run(path):
    """benchmark name (without the /iterations suffix) -> list of real_time of each repetition (us)."""
    with open(path) as f:
        d = json.load(f)
    out = {}
    unit_scale = {"ns": 1e-3, "us": 1.0, "ms": 1e3, "s": 1e6}
    for b in d["benchmarks"]:
        if b.get("run_type") != "iteration":
            continue
        name = re.sub(r"/iterations:\d+.*$", "", b["run_name"] if "run_name" in b else b["name"])
        out.setdefault(name, []).append(b["real_time"] * unit_scale[b.get("time_unit", "us")])
    return out, d.get("context", {})


def key_of(name):
    m = re.match(r"n=(\d+) / log\(q\)=(\d+) / (.*)$", name)
    n, op = int(m.group(1)), m.group(3)
    return n, op


def main():
    rawdir, outfile = sys.argv[1], sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) > 3 else ""
    mode_a = sys.argv[4] if len(sys.argv) > 4 else "off"
    mode_b = sys.argv[5] if len(sys.argv) > 5 else "on"
    pairs = {}  # (tag, group) -> {"off": run, "on": run}
    for fn in sorted(os.listdir(rawdir)):
        m = re.match(r"(.+)-(small|big)-(off|on|onscalar|offmain|offnear)\.json$", fn)
        if not m or not m.group(1).startswith(prefix):
            continue
        try:
            run, ctx = load_run(os.path.join(rawdir, fn))
        except (ValueError, KeyError) as e:  # an unfinished or broken file: say so, leave it out
            print(f"WARNING: {fn} left out: {e}", file=sys.stderr)
            continue
        if not run:
            continue
        pairs.setdefault((m.group(1), m.group(2)), {})[m.group(3)] = run
    complete = {k: v for k, v in pairs.items() if mode_a in v and mode_b in v}
    rows = {}  # (n, op) -> list of (off_med, on_med)
    reps = set()
    for (tag, group), pr in sorted(complete.items()):
        for name, offv in pr[mode_a].items():
            if name not in pr[mode_b]:
                continue
            onv = pr[mode_b][name]
            reps.add((len(offv), len(onv)))
            rows.setdefault(key_of(name), []).append((statistics.median(offv), statistics.median(onv), tag))
    npairs = {}
    for (tag, group) in complete:
        npairs[group] = npairs.get(group, 0) + 1
    lines = []
    lines.append(f"compared: {mode_a} (as OFF) with {mode_b} (as ON); pairs used: " + ", ".join(f"{g}: {c}" for g, c in sorted(npairs.items())) +
                 f"; repetitions per run (off, on): {sorted(reps)}")
    lines.append("| n | benchmark | OFF median (us) | ON median (us) | change | spread across pairs | verdict |")
    lines.append("|---|---|---|---|---|---|---|")
    csv = ["n,benchmark,pairs,off_median_us,on_median_us,change_pct,spread_lo_pct,spread_hi_pct,verdict"]
    for n in sorted({k[0] for k in rows}):
        for op in ORDER_OP + sorted({k[1] for k in rows if k[1] not in ORDER_OP}):
            if (n, op) not in rows:
                continue
            v = rows[(n, op)]
            offm = statistics.median(x[0] for x in v)
            onm = statistics.median(x[1] for x in v)
            ch = onm / offm - 1
            per = [x[1] / x[0] - 1 for x in v]
            lo, hi = min(per), max(per)
            if lo <= 0 <= hi:
                verdict = "no change"
            elif hi < 0:
                verdict = "faster"
            else:
                verdict = "slower"
            fmt = (lambda t: f"{t:,.1f}" if t < 1000 else f"{t:,.0f}")
            lines.append(f"| {n} | {op} | {fmt(offm)} | {fmt(onm)} | {ch*100:+.1f} % | {lo*100:+.1f} to {hi*100:+.1f} % "
                         f"| {verdict} |")
            csv.append(f"{n},{op},{len(v)},{offm:.4f},{onm:.4f},{ch*100:.4f},{lo*100:.4f},{hi*100:.4f},{verdict}")
    text = "\n".join(lines) + "\n"
    with open(outfile, "w") as f:
        f.write(text)
    with open(os.path.splitext(outfile)[0] + ".csv", "w") as f:
        f.write("\n".join(csv) + "\n")
    print(text, end="")


if __name__ == "__main__":
    main()
