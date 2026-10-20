"""coverage.py: own orbit canonicalisation (min 31-bit mask over all 155 maps k -> 2^j k + c), Burnside counts,
and coverage / verdict bookkeeping of the rm2-exact output files."""
import sys, glob, numpy as np
from math import comb
J = '../'
def burnside(W):
    fix = 0
    # a != 1: one fixed point + six 5-cycles
    for t in range(0, 7):
        for e in (0, 1):
            if 5 * t + e == W: fix += comb(6, t)
    return (comb(31, W) + 124 * fix) // 155
MAPS = [[(pow(2, j, 31) * k + c) % 31 for k in range(31)] for j in range(5) for c in range(31)]
def canon(masks):
    masks = np.asarray(masks, dtype=np.int64); best = np.full(len(masks), 1 << 40, dtype=np.int64)
    for mp in MAPS:
        img = np.zeros(len(masks), dtype=np.int64)
        for k in range(31): img |= ((masks >> k) & 1) << mp[k]
        best = np.minimum(best, img)
    return best
def read(files):
    rows = []
    for f in files:
        for l in open(f):
            t = l.split()
            if not t or t[0] != 'S': continue
            ks = [int(x) for x in t[1].split(',')]
            rows.append((sum(1 << k for k in ks), len(ks), len(set(ks)), t[2] if len(t) > 2 else '?', l.strip()))
    return rows
def report(name, rows, W):
    m = [r[0] for r in rows]; c = canon(m)
    assert all(r[1] == W and r[2] == W for r in rows)
    from collections import Counter
    print(name, "lines", len(rows), "distinct orbits", len(set(c.tolist())), "burnside", burnside(W), "verdicts", dict(Counter(r[3] for r in rows)))
    return dict(zip(c.tolist(), [r[3] for r in rows])), rows
if __name__ == "__main__":
    for W in (8, 9, 10): print("W", W, "burnside", burnside(W))
    d10, r10 = report("W10 s10", read(sorted(glob.glob(J + 'jobs/s10/remote-out/*/out_10_*.txt'))), 10)
    y10, _ = report("W10 s10y", read(sorted(glob.glob(J + 'jobs/s10y/remote-out/*/out_y_*.txt'))), 10)
    p10 = {k for k, v in d10.items() if v == 'P'}
    print("  s10 P == s10y set:", p10 == set(y10), "s10y all Y:", set(y10.values()))
    d9, _ = report("W9 s9y", read(sorted(glob.glob(J + 'jobs/s9y/remote-out/*/out_9_*.txt'))), 9)
    q9, _ = report("W9 q9", read(sorted(glob.glob(J + 'jobs/q9/remote-out/*/out_q9.txt'))), 9)
    p9 = {k for k, v in d9.items() if v == 'P'}
    print("  s9y P == q9 set:", p9 == set(q9), "q9 verdicts:", set(q9.values()))
    d8, r8 = report("W8 all8", read([J + 'all8.txt']), 8)
    q8lines = set()
    for f in glob.glob(J + 'jobs/q8*/remote-out/*/out_q8*.txt'):
        q8lines |= {l.strip() for l in open(f) if l.startswith('S ')}
    print("  all8 lines found verbatim in q8* outputs:", sum(r[4] in q8lines for r in r8), "of", len(r8))
    deep, _ = report("W8 deep*", read(sorted(glob.glob(J + 'jobs/deep*/remote-out/*/out_deep_*.txt'))), 8)
    p8 = {k for k, v in d8.items() if v == 'P'}
    print("  all8 P == deep set:", p8 == set(deep))
    fin = canon([sum(1 << int(x) for x in l.split(',')) for l in open(J + 'surv8_final.txt') if l.strip()]).tolist()
    pdeep = {k for k, v in deep.items() if v == 'P'}
    print("  deep P == surv8_final:", pdeep == set(fin), len(pdeep), len(fin))
    # mirror of 8a
    a = [0, 2, 4, 8, 18, 24, 25, 27]; mir = [(-k) % 31 for k in a]
    ca, cm = canon([sum(1 << k for k in a), sum(1 << k for k in mir)]).tolist()
    print("8a verdict", d8.get(ca), "mirror(-k) verdict", d8.get(cm), "mirror in final survivors:", cm in set(fin))
