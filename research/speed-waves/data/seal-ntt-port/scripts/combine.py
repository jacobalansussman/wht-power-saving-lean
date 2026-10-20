#!/usr/bin/env python3
"""combine.py : the tables of RESULTS.md. One table per polynomial size from m4-table.csv (7 OFF/ON pairs of the
main rounds), with the control experiment's same-binary figure (control-neon-alone.csv: the ON binary with the NEON
NTT running vs the same binary with SEAL's integer NTT running, 3 rounds) and a short note."""
import csv, os
H = os.path.dirname(os.path.abspath(__file__))
main = list(csv.DictReader(open(os.path.join(H, "m4-table.csv"))))
ctl = {(r["n"], r["benchmark"]): r for r in csv.DictReader(open(os.path.join(H, "control-neon-alone.csv")))}
NO_NTT = {"CKKS / EvaluateMulCt", "CKKS / EvaluateSquare", "CKKS / EvaluateAddCt", "CKKS / Decrypt"}
def note(n, b):
    if b == "CKKS / Decrypt":
        return "no NTT inside: code layout of the OFF binary, not the NEON code"
    if b in NO_NTT:
        return "no NTT inside (control)"
    if b == "BFV / Decrypt" and n != "32768":
        return "3-4 points of this are code layout"
    if n == "32768" and b in ("BFV / EncodeBatch", "BFV / DecodeBatch"):
        return "20-bit plaintext prime: NEON runs"
    if n == "32768":
        return "55-56-bit primes: NEON not used"
    return ""
def f(x):
    x = float(x)
    return f"{x:,.1f}" if x < 1000 else f"{x:,.0f}"
out = []
for n in ("4096", "8192", "16384", "32768"):
    out.append(f"\n#### n = {n}\n")
    out.append("| benchmark | OFF (µs) | ON (µs) | change | spread, 7 pairs | verdict | NEON code alone (control) | note |")
    out.append("|---|---|---|---|---|---|---|---|")
    for r in main:
        if r["n"] != n:
            continue
        c = ctl.get((n, r["benchmark"]))
        cs = f"{float(c['change_pct']):+.1f} % ({c['verdict']})" if c else "-"
        out.append(f"| {r['benchmark']} | {f(r['off_median_us'])} | {f(r['on_median_us'])} | {float(r['change_pct']):+.1f} % "
                   f"| {float(r['spread_lo_pct']):+.1f} to {float(r['spread_hi_pct']):+.1f} % | {r['verdict']} | {cs} | {note(n, r['benchmark'])} |")
open(os.path.join(H, "m4-tables-by-size.md"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
