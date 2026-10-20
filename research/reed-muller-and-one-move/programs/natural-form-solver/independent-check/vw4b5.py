"""exact tests (mine, incl. Z and TR) on a random sample of rm2-solve's block-5 W=4 S+D passes (their coordinates)."""
import sys, re, random, time
sys.path.insert(0, '..')
import packing2
from vtests import Tester
t = Tester(5)
lines = [l for l in open('../out_scan_frames4_partial.txt') if l.startswith('PASS')]
n = int(sys.argv[1]); seed = int(sys.argv[2])
random.Random(seed).shuffle(lines)
cnt = {"Z": 0, "TR": 0, "all": 0}; t0 = time.time()
for l in lines[:n]:
    mm = re.match(r"PASS ks=\(([^)]*)\) frames=\(([^)]*)\) delete dir (\d+)", l)
    ks = [int(v) for v in mm.group(1).split(',')]; ch = [int(v) for v in mm.group(2).split(',')]; z = int(mm.group(3))
    ms = [y for k, c in zip(ks, ch) for y in sorted(sorted(F) for F in packing2.all_frames(k))[c]]
    ms.remove(z)
    Ts = [packing2.fd.target(k) for k in ks]
    r = t.run(ms, None, Ts=Ts)
    assert r["S"] and r["D"], l
    cnt["all"] += 1; cnt["Z"] += r["C"] and r["Z"]; cnt["TR"] += r["C"] and r["Z"] and r["TR"]
    print(tuple(ks), tuple(ch), "del", z, {k: v for k, v in r.items() if k not in ("S", "D")}, f"{time.time()-t0:.0f}s", flush=True)
print("sampled", cnt["all"], "S+D passes; surviving C+Z:", cnt["Z"], "; C+Z+TR:", cnt["TR"])
