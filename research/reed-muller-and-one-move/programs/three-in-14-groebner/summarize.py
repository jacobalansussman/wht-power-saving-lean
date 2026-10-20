"""summarize.py: collect res_gb2_*.jsonl (fetched into jobs/gb2/remote-out/*/) into the final table."""
import json, glob, collections
rows = []
for f in sorted(glob.glob("jobs/gb2/remote-out/*/res_gb2_*.jsonl")):
    rows += [json.loads(l) for l in open(f)]
seen = {}
for r in rows: seen[r["idx"]] = r
rows = [seen[k] for k in sorted(seen)]
groups = collections.OrderedDict()
for r in rows:
    g = groups.setdefault(r["tag"], dict(n=0, orders=0, Qunit=0, Qnot=0, Qmiss=0, ms_unit=0, ms_not=0, ms_err=0,
                                         lifts_ok=0, lifts_checked=0, agree=0, expected=0, verdicts=collections.Counter(),
                                         dens=set(), secs=0.0))
    g["n"] += 1; g["orders"] += r["norders"]; g["Qunit"] += r["Q"]["unit"]; g["Qnot"] += r["Q"]["not_unit"]
    g["Qmiss"] += r["Q"]["missing"]
    for p, v in r["msolve"].items():
        g["ms_unit"] += v["unit"]; g["ms_not"] += v["not_unit"]; g["ms_err"] += v["error"]
    g["lifts_ok"] += r["lifts"].get("ok", 0); g["lifts_checked"] += r["lifts"].get("checked", 0)
    g["dens"] |= set(r["lifts"].get("denominators", []))
    g["agree"] += r["agree"]; g["expected"] += r["as_expected"]; g["verdicts"][r["verdict"]] += 1; g["secs"] += r["secs"]
print("| group | instances | orders | Singular/Q: 1 in ideal / not | msolve (2 primes): unit / not / error | lift identities checked exactly | all three agree | as expected | verdicts |")
print("|---|---|---|---|---|---|---|---|---|")
for t, g in groups.items():
    print(f"| {t} | {g['n']} | {g['orders']} | {g['Qunit']} / {g['Qnot']}" + (f" (missing {g['Qmiss']})" if g['Qmiss'] else "") +
          f" | {g['ms_unit']} / {g['ms_not']} / {g['ms_err']} | {g['lifts_ok']} of {g['lifts_checked']} | {g['agree']} of {g['n']} | "
          f"{g['expected']} of {g['n']} | {dict(g['verdicts'])} |")
print()
for t, g in groups.items():
    print(t, "lift denominators:", sorted(g["dens"])[:20], "cpu-s", round(g["secs"]))
bad = [r for r in rows if not r["as_expected"]]
print("instances not as expected:", [(r["idx"], r["tag"], r["verdict"]) for r in bad])
print("total instances", len(rows))
