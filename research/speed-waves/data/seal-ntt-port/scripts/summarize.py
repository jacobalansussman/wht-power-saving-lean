#!/usr/bin/env python3
"""summarize.py [--csv FILE] [--compact] JOB_DIR [JOB_DIR ...]

Reads what one or more jobs left in their output folders (sysinfo.txt, tests-*.status, ab-meta.json and the
Google Benchmark files bench-{off,on}-pN.json written by ab.sh) and prints Markdown tables.

For each benchmark and each pair, the time of a run is the median of its repetitions (wall-clock time, in
microseconds). Over the pairs, the table gives the median OFF time, the median ON time, the change of the medians
(ON / OFF - 1), and the spread: the smallest and largest change seen within a single pair.
Standard library only.
"""
import csv
import glob
import json
import os
import re
import statistics
import sys

OPS = {
    ("UTIL", "NTTForward"): "NTT forward, all primes",
    ("UTIL", "NTTInverse"): "NTT inverse, all primes",
    ("UTIL", "NTTForwardLowLevel"): "NTT forward, one prime",
    ("UTIL", "NTTInverseLowLevel"): "NTT inverse, one prime",
    ("UTIL", "NTTForwardLowLevelLazy"): "NTT forward, one prime, lazy",
    ("UTIL", "NTTInverseLowLevelLazy"): "NTT inverse, one prime, lazy",
    ("BFV", "EncryptSecret"): "BFV encrypt (secret key)",
    ("BFV", "EncryptPublic"): "BFV encrypt (public key)",
    ("BFV", "Decrypt"): "BFV decrypt",
    ("BFV", "EvaluateMulCt"): "BFV multiply",
    ("BFV", "EvaluateRelinInplace"): "BFV relinearize",
    ("BFV", "EvaluateRotateRows"): "BFV rotate rows",
    ("CKKS", "EncryptSecret"): "CKKS encrypt (secret key)",
    ("CKKS", "EncryptPublic"): "CKKS encrypt (public key)",
    ("CKKS", "Decrypt"): "CKKS decrypt",
    ("CKKS", "EvaluateMulCt"): "CKKS multiply (no NTT inside)",
    ("CKKS", "EvaluateRescaleInplace"): "CKKS rescale",
    ("CKKS", "EvaluateRelinInplace"): "CKKS relinearize",
    ("CKKS", "EvaluateRotate"): "CKKS rotate",
}
ORDER = list(OPS.values())
NAME = re.compile(r"^n=(\d+) / log\(q\)=(\d+) / (\w+) / (\w+)/")


def parse_name(run_name):
    m = NAME.match(run_name)
    if not m:
        return None
    n, _logq, scheme, op = m.groups()
    return int(n), OPS.get((scheme, op), scheme + " " + op)


def load_runs(job, mode):
    """{pair: {(n, op): median real_time in us}} for one build."""
    out = {}
    for f in glob.glob(os.path.join(job, "bench-%s-p*.json" % mode)):
        pair = int(re.search(r"-p(\d+)\.json$", f).group(1))
        try:
            with open(f) as fh:
                data = json.load(fh)
        except (OSError, ValueError):
            continue
        times = {}
        for b in data.get("benchmarks", []):
            if b.get("run_type") == "aggregate" and b.get("aggregate_name") == "median":
                key = parse_name(b["run_name"])
                if key:
                    unit = b.get("time_unit", "us")
                    scale = {"ns": 1e-3, "us": 1.0, "ms": 1e3, "s": 1e6}[unit]
                    times[key] = b["real_time"] * scale
        out[pair] = times
    return out


def job_rows(job):
    off, on = load_runs(job, "off"), load_runs(job, "on")
    pairs = sorted(set(off) & set(on))
    keys = set()
    for p in pairs:
        keys |= set(off[p]) & set(on[p])
    rows = []
    for key in sorted(keys, key=lambda k: (ORDER.index(k[1]) if k[1] in ORDER else 99, k[0])):
        o = [off[p][key] for p in pairs if key in off[p] and key in on[p]]
        a = [on[p][key] for p in pairs if key in off[p] and key in on[p]]
        if not o:
            continue
        per_pair = [100.0 * (y / x - 1.0) for x, y in zip(o, a)]
        mo, ma = statistics.median(o), statistics.median(a)
        rows.append({
            "n": key[0], "operation": key[1], "pairs": len(o),
            "off_us": mo, "on_us": ma, "change_pct": 100.0 * (ma / mo - 1.0),
            "pair_min_pct": min(per_pair), "pair_max_pct": max(per_pair),
        })
    return rows, len(pairs)


def fmt_us(x):
    if x >= 1000:
        return "{:,.0f}".format(x)
    if x >= 100:
        return "{:.0f}".format(x)
    if x >= 10:
        return "{:.1f}".format(x)
    return "{:.2f}".format(x)


def fmt_pct(x):
    return ("+" if x >= 0.05 else ("−" if x <= -0.05 else "")) + "{:.1f} %".format(abs(x))


def machine(job):
    info = {}
    try:
        with open(os.path.join(job, "sysinfo.txt")) as fh:
            text = fh.read()
    except OSError:
        return info
    m = re.search(r"^Model name:\s*(.+)$", text, re.M)
    info["cpu"] = m.group(1).strip() if m else "?"
    m = re.search(r"^CPU\(s\):\s*(\d+)$", text, re.M)
    info["cpus"] = m.group(1) if m else "?"
    m = re.search(r"^### \S+ --version\n(.+)$", text, re.M)
    info["cc"] = m.group(1).strip() if m else "?"
    m = re.search(r"^MemTotal:\s*(\d+) kB", text, re.M)
    info["mem_gb"] = "%.1f" % (int(m.group(1)) / 1048576) if m else "?"
    return info


def tests(job):
    lines = []
    for f in sorted(glob.glob(os.path.join(job, "tests-*.status"))):
        with open(f) as fh:
            lines.append(dict(kv.split("=", 1) for kv in fh.read().split()))
    return lines


def print_job(job):
    name = os.path.basename(os.path.normpath(job))
    mi = machine(job)
    print("### %s\n" % name)
    if mi:
        print("Machine: %s, %s CPUs, %s GB; compiler: %s\n" % (mi["cpu"], mi["cpus"], mi["mem_gb"], mi["cc"]))
    t = tests(job)
    if t:
        print("| build | tests passed | failed | exit code | seconds |")
        print("|---|---|---|---|---|")
        for s in t:
            print("| %s | %s of %s | %s | %s | %s |" % (s["mode"], s["passed"], s["listed"], s["failed"], s["exit"],
                                                      s["seconds"]))
        print()
    rows, npairs = job_rows(job)
    if rows:
        print("%d pairs. Times in microseconds per operation (median over pairs of each run's median)." % npairs)
        print()
        print("| operation | n | OFF | ON | change | spread over pairs |")
        print("|---|---:|---:|---:|---:|---|")
        for r in rows:
            print("| %s | %d | %s | %s | %s | %s to %s |" % (
                r["operation"], r["n"], fmt_us(r["off_us"]), fmt_us(r["on_us"]), fmt_pct(r["change_pct"]),
                fmt_pct(r["pair_min_pct"]), fmt_pct(r["pair_max_pct"])))
        print()
    return rows


def print_compact(jobs):
    """One table: change per operation and size, one column group per job."""
    data = {}
    for job in jobs:
        rows, _ = job_rows(job)
        if rows:
            data[os.path.basename(os.path.normpath(job))] = {(r["operation"], r["n"]): r for r in rows}
    if not data:
        return
    names = list(data)  # in the order given on the command line
    sizes = sorted({k[1] for d in data.values() for k in d})
    ops = [o for o in ORDER if any(o == k[0] for d in data.values() for k in d)]
    print("| operation | " + " | ".join("n = %d" % n for n in sizes) + " |")
    print("|---|" + "---:|" * len(sizes))
    for op in ops:
        cells = []
        for n in sizes:
            vals = [fmt_pct(data[j][(op, n)]["change_pct"]) if (op, n) in data[j] else "–" for j in names]
            cells.append(" / ".join(vals))
        print("| %s | %s |" % (op, " | ".join(cells)))
    print("\nEach cell: change of the ON build against the OFF build, " + " / ".join(names) + ".\n")


def main(argv):
    csv_path, compact, jobs = None, False, []
    it = iter(argv)
    for a in it:
        if a == "--csv":
            csv_path = next(it)
        elif a == "--compact":
            compact = True
        else:
            jobs.append(a)
    if not jobs:
        print(__doc__)
        return 2
    if compact:
        print_compact(jobs)
    all_rows = []
    for job in jobs:
        for r in print_job(job):
            r["job"] = os.path.basename(os.path.normpath(job))
            all_rows.append(r)
    if csv_path and all_rows:
        with open(csv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["job", "operation", "n", "pairs", "off_us", "on_us", "change_pct",
                                               "pair_min_pct", "pair_max_pct"])
            w.writeheader()
            w.writerows(all_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
