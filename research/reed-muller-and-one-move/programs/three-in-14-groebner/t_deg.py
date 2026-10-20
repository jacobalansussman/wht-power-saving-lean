import json, subprocess, os, sys
from f2q import Field, Config
from gb2 import supports_in, orders, Ring, system, pstr
certs = json.load(open("certs_min.json")); F = Field(5)
for idx in (0, 2, 88, 132):
    c = certs[idx]; cfg = Config(Field(c["m"]), c["Ts"], c["moves"])
    U = [i for i in range(cfg.M) if c["cert"]["U"] >> i & 1]; sup = supports_in(cfg, U); R = Ring(cfg.W, U)
    o = list(orders(cfg, U))[0]; ps = system(cfg, R, o, sup); NV = len(R.names)
    for d in (4, 5, 6, 7, 8, 10):
        txt = f"ring r = 32003, (x(1..{NV})), dp;\ndegBound = {d};\nideal I = " + ",\n".join(pstr(p, "sing") for p in ps) + ";\nideal G = std(I); if (G[1]==1) {\"UNIT\";} else {\"NF\";}\nquit;\n"
        open("tdeg.sing", "w").write(txt)
        out = subprocess.run(["Singular", "-q", "tdeg.sing"], capture_output=True, text=True, timeout=120).stdout.strip()
        print(idx, c["tag"], "degBound", d, out, flush=True)
        if out == "UNIT": break
