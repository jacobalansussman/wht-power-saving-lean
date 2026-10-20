"""z4b5.py: my span/deletion/clique/zero-block tests on v-rm2-solve's nine block-5 4-in-19 survivors."""
import re
from f2q import Field, Config
from ztest import zero_block
F = Field(5)
for k, line in enumerate(l for l in open("survivors_b5w4.txt") if l.startswith("SURVIVES")):
    ks = [int(v) for v in re.search(r"\(([0-9, ]+)\)", line).group(1).split(",")]
    moves = [int(v) for v in re.search(r"dirs \[([0-9, ]+)\]", line).group(1).split(",")]
    cfg = Config(F, [F.target(x) for x in ks], moves); r = zero_block(cfg)
    print(k, ks, "span/deletion", cfg.span_deletion(), "clique", r["clique"], "C", r["C"], "Z", r["Z"], "tight supports", r["ntight_supports"], flush=True)
