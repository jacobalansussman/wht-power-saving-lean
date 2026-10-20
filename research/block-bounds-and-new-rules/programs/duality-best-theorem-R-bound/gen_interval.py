# gen_interval.py m -> interval_m{m}.txt : frames N_s = s (+) s^0 (basis vectors as 2m-bit ints) and the
# unordered edges / complementary pairs of the interval (pf-duality).  Read by anneal_iv.c.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lift_interval import frames, edges
from interval import dimof
m = int(sys.argv[1])
fr = frames(m); subs = [s for s, _ in fr]
E = edges(subs)
comp = [(i, j) for i, s in enumerate(subs) for j, t in enumerate(subs)
        if i < j and len(s & t) == 1 and dimof(s) + dimof(t) == m]
with open(f"interval_m{m}.txt", "w") as f:
    f.write(f"{m} {len(fr)} {len(E)} {len(comp)}\n")
    for s, vecs in fr: f.write(" ".join(map(str, vecs)) + "\n")
    for i, j in E: f.write(f"{i} {j}\n")
    for i, j in comp: f.write(f"{i} {j}\n")
print(m, len(fr), "frames", len(E), "edges", len(comp), "complementary pairs")
