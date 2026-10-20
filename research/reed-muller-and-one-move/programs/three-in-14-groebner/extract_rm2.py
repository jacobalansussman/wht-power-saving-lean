"""extract_rm2.py: DATA EXTRACTION ONLY. Uses rm2-solve's own field/frames (its coordinates) to turn each PASS line
of rm2-solve/out_scan_frames3.txt into (targets as 0/1 matrices, the 14 moves as direction ints), written to
data_rm2_44.json. All tests on these instances are done by my own code (ztest.py, gb.py)."""
import sys, re, json, os
sys.path.insert(0, "../natural-form-solver")
os.environ["HM"] = "5"
from packing2 import all_frames, fd
out = []
for line in open("../natural-form-solver/out_scan_frames3.txt"):
    if not line.startswith("PASS"): continue
    m = re.search(r"ks=\(([0-9, ]+)\) frames=\(([0-9, ]+)\) delete dir (\d+)", line)
    ks = [int(v) for v in m.group(1).split(",")]; ch = [int(v) for v in m.group(2).split(",")]; z = int(m.group(3))
    plain = [w for k, c in zip(ks, ch) for w in sorted(sorted(F) for F in all_frames(k))[c]]
    d = list(plain); d.remove(z)
    out.append(dict(line=line.strip(), ks=ks, frames=ch, delete=z, Ts=[fd.target(k) for k in ks], plain=plain, moves=d))
json.dump(out, open("data_rm2_44.json", "w"))
print(len(out), "instances written")
