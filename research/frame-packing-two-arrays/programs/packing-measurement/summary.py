#!/usr/bin/env python3
"""fp-measure: one table of the raw-gauge ledger of every object measured (reads out/*.json written by calib.py,
fm_ledger.py, oa_ledger.py, oa_closed.py), plus the whole-unit estimates for the units whose global word was not
replayed (E8 unit, U8X0HQ, T74, the block-110 unit): calls x (measured forward call) / (the unit's moves).
usage: python3 -I summary.py > out/summary.txt"""
import json, os, math

O = "out"
rows = []
c = json.load(open(os.path.join(O, "calib.json")))
seen = set()
for r in c["a"]:
    k = (r["case"], r["m"], r["W"])
    if k in seen:
        continue
    seen.add(k)
    rows.append(dict(obj=r["case"], kind="unit, exact (engine)", m=r["m"], W=r["W"], M=r["M"], Delta=r["Delta"], P=r["P"],
                     B=r["B"], D=r["D"], Z=r["Z"], I0=r["I0"], maxlab=r["maxlab"]))
for r in json.load(open(os.path.join(O, "oa_ledger_6_7_8_9.json"))):
    rows.append(dict(obj="OpenAI network h=%d" % r["h"], kind="unit, exact (frame model)", m=r["m"], W=r["W"], M=r["M"],
                     Delta=r["Delta"], P=r["P"], B=r["B"], D=r["D"], Z=0, I0=0, maxlab=r["maxlab"]))
for r in json.load(open(os.path.join(O, "oa_closed.json"))):
    if r["h"] in (22, 25, 100):
        rows.append(dict(obj="OpenAI network h=%d" % r["h"], kind="unit, closed form", m=r["m"], W=r["W"], M=r["M"],
                         Delta=r["Delta"], P=r["P"], B=r["B"], D=r["D"], Z=0, I0=0, maxlab=r["maxlab"]))
calls = {"E8r783": ("E8 unit invocation (h=9, v=120)", 45), "U8X": ("U8 line call (U8X, h=9, v=112)", 73),
         "X4X": ("#799 call (X4X, h=9, v=112)", 73), "T74g": ("T74 line call (gauge44-566, h=9, v=112)", 73),
         "P11pr193": ("block-110 invocation (pr193, h=22, v=1320)", 110)}
per_call = {}
for n, (label, mm) in calls.items():
    p = os.path.join(O, "fm_%s.json" % n)
    if not os.path.exists(p):
        continue
    d = json.load(open(p))
    r = d["modes"]["none"]
    per_call[n] = d["modes"]
    for mode, tag in (("none", "F1 fresh"), ("F2/none", "F2 src x-y"), ("F3/none", "F3 src -y")):
        if mode not in d["modes"]:
            continue
        r = d["modes"][mode]
        rows.append(dict(obj=label, kind="one call %s, exact" % tag, m=mm, W=2 * d["v"] + d["R"], M=r["M"],
                         Delta=r["deliv"] - r["M"], P=r["P"], B=r["B"], D=r["D"], Z=r["Z"], I0=r["I0"],
                         maxlab=r["maxlab"]))
print("%-44s %-30s %8s %10s %12s %13s %14s %14s %13s %6s %4s | %8s %7s %8s %8s" % (
    "object", "what", "block", "arrays", "moves M", "Delta", "P", "B", "D", "I0", "max",
    "P/M", "B/P", "2D/M", "(P-B)/M"))
for r in rows:
    M = r["M"]
    print("%-44s %-30s %8d %10.4g %12.6g %13.6g %14.6g %14.6g %13.6g %6d %4d | %8.4f %7.4f %8.4f %8.4f" % (
        r["obj"], r["kind"], r["m"], r["W"], M, r["Delta"], r["P"], r["B"], r["D"], r["I0"], r["maxlab"],
        r["P"] / M, (r["B"] / r["P"]) if r["P"] else float("nan"), 2 * r["D"] / M, (r["P"] - r["B"]) / M))
print("\n(For one call, 'Delta' is the call's own balance: delivered levels - moves; the unit's Delta is not a sum of these.)")
# per-call matrix: P/M by kind of call and helper model
print("\nP per move by kind of call (F1 fresh forward, F2 forward with sources x - y, F3 forward with sources -y, B1 backward"
      " with fresh sources, B2 backward with sources y + x) and helper model (clean | dirty, one dirt class | dirty, a class"
      " per helper):")
FW = {"clean": "", "one": "/gates-one", "each": "/gates"}
BW = {"clean": "/clean", "one": "/dirty-one", "each": "/dirty"}
callP = {}
for n in ("E8r783", "U8X", "X4X", "T74g", "P11pr193"):
    pf, pb = os.path.join(O, "fm_%s.json" % n), os.path.join(O, "bw_%s.json" % n)
    if not os.path.exists(pf):
        continue
    f = json.load(open(pf))["modes"]
    b = json.load(open(pb)) if os.path.exists(pb) else {}
    callP[n] = {}
    print("  %s" % calls[n][0])
    for kind in ("F1", "F2", "F3", "B1", "B2"):
        cells = []
        for hm in ("clean", "one", "each"):
            if kind[0] == "F":
                key = ("none" if hm == "clean" else FW[hm][1:]) if kind == "F1" else kind + ("/none" if hm == "clean" else FW[hm])
                r = f.get(key)
            else:
                r = b.get(kind + BW[hm])
            if r:
                callP[n][(kind, hm)] = (r["P"], r["M"])
                cells.append("%7.2f" % (r["P"] / r["M"]))
            else:
                cells.append("    n/a")
        print("    %s  %s" % (kind, " | ".join(cells)))
est = [("E8 unit (m=45)", "E8r783", 5, 1263 * 45 - 120, 120),
       ("U8X0HQ (m=73, 8 replicas)", "U8X", 72, 73 * 7491 - 2560, 2560),
       ("T74 per replica (m=73)", "T74g", 9, 73 * 896 - 320, 320),
       ("block-110 unit (m=110)", "P11pr193", 5, 14692 * 110 - 3080, 3080)]
print("\nWhole-unit ESTIMATES. Global words not replayed (gluing, idle climbs, hosting, re-framing). Every call taken as one"
      " shear of a three-shear swap, F1, B2, F3 in equal shares (B2 replaced by F1 where the backward call was not run):")
for name, n, k, M, Dl in est:
    if n not in callP:
        continue
    out = []
    for hm in ("clean", "one", "each"):
        c = callP[n]
        if ("F1", hm) not in c or ("F3", hm) not in c:
            out.append("n/a"); continue
        b2 = c.get(("B2", hm), c[("F1", hm)])
        avg = (c[("F1", hm)][0] + b2[0] + c[("F3", hm)][0]) / 3
        out.append("%.2f" % (k * avg / M))
    print("  %-28s calls %3d  M %9d  P/M ~ %s (clean | dirty-one | dirty-each)   Delta/M = %.5f  (Delta = %d)" % (
        name, k, M, " | ".join(out), Dl / M, Dl))
