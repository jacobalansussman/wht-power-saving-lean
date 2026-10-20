"""Tables for RESULTS.md: per compiler, run A in full plus run B's change; a compact table; an M4 comparison.
usage: python3 combine_runs.py SCRIPTS_DIR RUN_A_DIR RUN_B_DIR M4_DIR
  (folders of raw benchmark files per job; the raw files are not in this revision)"""
import sys

sys.path.insert(0, sys.argv[1])
import summarize as S  # noqa: E402

A, B, M4 = sys.argv[2], sys.argv[3], sys.argv[4]
comp = [("GCC 13.3", "arm64-gcc"), ("Clang 18.1", "arm64-clang")]


def rows(d):
    r, n = S.job_rows(d)
    return {(x["operation"], x["n"]): x for x in r}, n


ra = {j: rows(A + "/" + j) for _, j in comp}
rb = {j: rows(B + "/" + j) for _, j in comp}
m4, _ = rows(M4)
ops = S.ORDER
sizes = [4096, 8192, 16384, 32768]

print("## COMPACT (run A / run B per compiler)\n")
print("| operation | n | GCC | GCC repeat | Clang | Clang repeat |")
print("|---|---:|---:|---:|---:|---:|")
for op in ops:
    for n in sizes:
        k = (op, n)
        cells = []
        for _, j in comp:
            for src in (ra, rb):
                cells.append(S.fmt_pct(src[j][0][k]["change_pct"]) if k in src[j][0] else "–")
        print("| %s | %d | %s |" % (op, n, " | ".join(cells)))
print()

for name, j in comp:
    print("## FULL %s (run A %d pairs, run B %d pairs)\n" % (name, ra[j][1], rb[j][1]))
    print("| operation | n | OFF (µs) | ON (µs) | change | spread over pairs | change, repeat run |")
    print("|---|---:|---:|---:|---:|---|---:|")
    for op in ops:
        for n in sizes:
            k = (op, n)
            if k not in ra[j][0]:
                continue
            x = ra[j][0][k]
            rep = S.fmt_pct(rb[j][0][k]["change_pct"]) if k in rb[j][0] else "–"
            print("| %s | %d | %s | %s | %s | %s to %s | %s |" % (
                op, n, S.fmt_us(x["off_us"]), S.fmt_us(x["on_us"]), S.fmt_pct(x["change_pct"]),
                S.fmt_pct(x["pair_min_pct"]), S.fmt_pct(x["pair_max_pct"]), rep))
    print()

print("## M4 vs N2 (run A)\n")
pick = ["NTT forward, one prime", "NTT inverse, one prime", "NTT forward, one prime, lazy",
        "NTT inverse, one prime, lazy", "BFV decrypt", "BFV multiply", "CKKS encrypt (public key)", "CKKS rescale",
        "CKKS relinearize", "CKKS rotate"]
print("| operation | n | Apple M4 (Apple clang) | Neoverse N2, GCC | Neoverse N2, Clang |")
print("|---|---:|---:|---:|---:|")
for op in pick:
    for n in [4096, 8192, 16384]:
        k = (op, n)
        a = S.fmt_pct(m4[k]["change_pct"]) if k in m4 else "–"
        g = S.fmt_pct(ra["arm64-gcc"][0][k]["change_pct"])
        c = S.fmt_pct(ra["arm64-clang"][0][k]["change_pct"])
        print("| %s | %d | %s | %s | %s |" % (op, n, a, g, c))
print()
print("## ON times GCC vs Clang, and OFF times (run A)\n")
for op in ["NTT forward, one prime", "NTT inverse, one prime", "CKKS relinearize", "BFV relinearize"]:
    for n in [4096, 16384]:
        k = (op, n)
        g, c = ra["arm64-gcc"][0][k], ra["arm64-clang"][0][k]
        print(op, n, "OFF gcc %.1f clang %.1f | ON gcc %.1f clang %.1f" % (g["off_us"], c["off_us"], g["on_us"],
                                                                        c["on_us"]))

print("## HEADLINE (run A; n = 4096, 8192, 16384; control n = 32768)\n")
head = ["NTT forward, one prime", "NTT inverse, one prime", "NTT forward, one prime, lazy",
        "NTT inverse, one prime, lazy", "BFV encrypt (public key)", "BFV decrypt", "BFV multiply", "BFV relinearize",
        "BFV rotate rows", "CKKS encrypt (public key)", "CKKS rescale", "CKKS relinearize", "CKKS rotate",
        "CKKS multiply (no NTT inside)"]


def trio(src, op):
    vals = []
    for n in [4096, 8192, 16384]:
        k = (op, n)
        if k not in src:
            return "–"
        vals.append(S.fmt_pct(src[k]["change_pct"]).replace(" %", ""))
    return ", ".join(vals) + " %"


print("| operation | GCC | Clang | Apple M4 (step 1) | n = 32768, GCC / Clang |")
print("|---|---:|---:|---:|---:|")
for op in head:
    g, c = ra["arm64-gcc"][0], ra["arm64-clang"][0]
    ctl = "%s / %s" % (S.fmt_pct(g[(op, 32768)]["change_pct"]), S.fmt_pct(c[(op, 32768)]["change_pct"]))
    print("| %s | %s | %s | %s | %s |" % (op, trio(g, op), trio(c, op), trio(m4, op), ctl))
print()
