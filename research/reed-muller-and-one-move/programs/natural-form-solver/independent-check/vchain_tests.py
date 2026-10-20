"""exact tests (mine) on the 8-in-39 deletion candidates of the rm2-solve chains, in rm2-solve's coordinates."""
import sys
import numpy as np
sys.path.insert(0, '..')
from inst import canonical
from vtests import Tester
t = Tester(5)
for spec in ["0,2,4,8,18,24,25,27", "0,0,1,9,16,18,21,25"]:
    ks = [int(v) for v in spec.split(',')]
    inst, fd = canonical(5, ks, dirs="frames", order="rand1", seed=0)
    Ts = inst.T; dirs = list(inst.dirs)
    base = t.run(dirs, None, Ts=Ts)
    print(spec, "plain 8 in 40:", {k: v for k, v in base.items()}, "dirs", dirs)
    seen = {}
    for i, z in enumerate(dirs):
        if z in seen: continue
        d2 = dirs[:i] + dirs[i + 1:]
        r = t.run(d2, None, Ts=Ts); seen[z] = r
        if r["S"] and r["D"]:
            print("  delete dir", z, "(mult", dirs.count(z), "):", {k: v for k, v in r.items() if k not in ("S", "D")}, flush=True)
    print("  directions whose deletion passes S+D:", [z for z, r in seen.items() if r["S"] and r["D"]])
