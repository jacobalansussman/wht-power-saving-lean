# summarise the exhaustive sweeps (out_box2_b4.txt, out_box2_b5_sep.txt): per (m, W, M) the number of deletion-passing
# multisets, how many the clique test excludes, survivors, and the saving the target would have.
import re, math, sys
rows = []
for fn in sys.argv[1:]:
    for l in open(fn):
        mm = re.match(r"m=(\d+) W=(\d+) M=(\d+): .*deletion-passing multisets (\d+); EXCLUDED by clique (\d+); survivors (\d+)", l)
        if mm:
            m, W, M, npass, nex, ns = map(int, mm.groups()); rows.append((m, W, M, npass, nex, ns))
        mm = re.match(r"m=(\d+) W=(\d+) M=(\d+): (LP infeasible|empty LP box)", l)
        if mm:
            m, W, M = map(int, mm.groups()[:3]); rows.append((m, W, M, 0, 0, 0))
print("m  W   M  saving  deletion-passers  excluded  survivors")
for m, W, M, a, b, c in sorted(set(rows)):
    if M >= m * W: continue
    print(f"{m}  {W:2d} {M:3d}  {math.log(m*W/M)/math.log(m):.3f}   {a:6d}            {b:6d}    {c:4d}")
