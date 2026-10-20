import sys, re
sys.path.insert(0, '..')
import packing2
from vtests import Tester
t = Tester(5)
for l in open('../cases4.txt'):
    mm = re.match(r"PASS ks=\(([^)]*)\) frames=\(([^)]*)\) delete dir (\d+)", l)
    if not mm: continue
    ks = [int(v) for v in mm.group(1).split(',')]; ch = [int(v) for v in mm.group(2).split(',')]; z = int(mm.group(3))
    ms = [y for k, c in zip(ks, ch) for y in sorted(sorted(F) for F in packing2.all_frames(k))[c]]
    ms.remove(z)
    r = t.run(ms, None, Ts=[packing2.fd.target(k) for k in ks])
    print(tuple(ks), tuple(ch), "del", z, {k: v for k, v in r.items()}, flush=True)
