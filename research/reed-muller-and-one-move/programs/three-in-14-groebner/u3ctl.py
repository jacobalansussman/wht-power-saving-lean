"""u3ctl.py: the must-pass control U_3 = 2 in 7 (block 3, targets I, I, all seven directions, lex order 1..7)
through the same exact algebra as the 44 (gb2.system with U = all 7 moves): Singular std over Q and msolve over
65521. Expected: NOT-UNIT (rm2-solve's solver finds U_3; ex2-families found it in lex order)."""
import subprocess, os, time
from f2q import Field, Config
from gb2 import supports_in, Ring, system, singular_text, msolve_text
F = Field(3); I3 = [[int(i == j) for j in range(3)] for i in range(3)]
cfg = Config(F, [I3, I3], list(range(1, 8))); U = list(range(7))
sup = supports_in(cfg, U); R = Ring(cfg.W, U); ps = system(cfg, R, U, sup)
print("U_3 lex: supports", len(sup), "polys", len(ps), "unknowns used", len({v for p in ps for m in p for v in m}), flush=True)
open("u3_Q.sing", "w").write(singular_text(R, [ps], 0, "U3lex"))
open("u3.ms", "w").write(msolve_text(ps, 65521))
t = time.time()
try:
    out = subprocess.run(["msolve", "-f", "u3.ms", "-o", "u3.out", "-g", "1"], capture_output=True, timeout=1500)
    s = open("u3.out").read(); print("msolve 65521:", "UNIT" if "[1]:" in s.split("#---")[-1][:20] else s.split("#length of basis:")[-1].split("\n")[0], f"{time.time() - t:.0f}s", flush=True)
except subprocess.TimeoutExpired:
    print("msolve timeout", flush=True)
t = time.time()
try:
    out = subprocess.run(["Singular", "-q", "u3_Q.sing"], capture_output=True, text=True, timeout=1800).stdout
    print("Singular Q:", [l for l in out.split("\n") if "ORDER" in l or "SUMMARY" in l], f"{time.time() - t:.0f}s", flush=True)
except subprocess.TimeoutExpired:
    print("Singular timeout", flush=True)
