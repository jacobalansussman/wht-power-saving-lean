"""z4b4.py: my zero-block test on my own block-4 4-in-15 near-plain list (out_scan4_b4.txt)."""
import re
from f2q import Field, Config, popc
from ztest import zero_block, check_cert
F = Field(4); n = 0; surv = []
for line in open("out_scan4_b4.txt"):
    if not line.startswith("PASS"): continue
    ks = [int(v) for v in re.search(r"ks=\(([0-9, ]+)\)", line).group(1).split(",")]
    moves = [int(v) for v in re.search(r"moves=\[([0-9, ]+)\]", line).group(1).split(",")]
    cfg = Config(F, [F.target(k) for k in ks], moves); r = zero_block(cfg); n += 1
    if r.get("cert") and r["cert"]["kind"] == "Z": check_cert(cfg, r["cert"])
    if r["C"] and r["Z"]: surv.append(line.strip()); print("SURVIVES C+Z", line.strip(), flush=True)
print(f"block 4, W=4: {n} passing span+deletion, {len(surv)} survive C+Z")
