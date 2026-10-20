"""summarize.py -- o3-shapeR: aggregate every family run (out/*.jsonl and remote-out/*/out/*.jsonl, deduplicated by file
name) into one table: prefixes, target pairs, PASS pairs, max dim eps(M4) / dim eps(M5') per prefix, best relaxed level."""
import json, glob, os, collections
H = os.path.dirname(os.path.abspath(__file__))
files = {}
for f in glob.glob(os.path.join(H, "out", "*.jsonl")) + glob.glob(os.path.join(H, "remote-out", "*", "out", "*.jsonl")):
    b = os.path.basename(f)
    if not b.startswith(("m1_", "m2_", "w1_", "w2_", "w3_", "w4_")):
        continue
    n = sum(1 for _ in open(f))
    if b not in files or n > files[b][1]:
        files[b] = (f, n)
tab = collections.defaultdict(lambda: collections.Counter())
for b, (f, n) in sorted(files.items()):
    for l in open(f):
        r = json.loads(l)
        fam = r["fam"]; T = tab[fam]
        T["prefixes"] += 1; T["pairs"] += r["npairs"]; T["pass pairs"] += r["npass"]
        T[("d4max", max((d[0] for d in r["d4"]), default=-1))] += 1
        T[("d5max", max((d[0] for d in r["d5"]), default=-1))] += 1
        T[("bestlev", r.get("bestlev", "n/a"))] += 1
print("files:", ", ".join(f"{b} ({n})" for b, (f, n) in sorted(files.items())))
tot = collections.Counter()
for fam in sorted(tab):
    T = tab[fam]
    d4 = {k[1]: v for k, v in T.items() if isinstance(k, tuple) and k[0] == "d4max"}
    d5 = {k[1]: v for k, v in T.items() if isinstance(k, tuple) and k[0] == "d5max"}
    bl = {str(k[1]): v for k, v in T.items() if isinstance(k, tuple) and k[0] == "bestlev"}
    print(f"{fam:7s} prefixes {T['prefixes']:5d} pairs {T['pairs']:6d} PASS {T['pass pairs']:3d} | max dim eps(M4) {dict(sorted(d4.items()))} | max dim eps(M5') {dict(sorted(d5.items()))} | best relaxed level {dict(sorted(bl.items()))}")
    tot["prefixes"] += T["prefixes"]; tot["pairs"] += T["pairs"]; tot["pass"] += T["pass pairs"]
print("TOTAL prefixes", tot["prefixes"], "target pairs", tot["pairs"], "PASS pairs", tot["pass"])
