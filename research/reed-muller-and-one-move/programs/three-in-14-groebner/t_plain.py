import json, re, os
from f2q import Field, Config
from gb2 import supports_in, orders, Ring, system, singular_text, msolve_text
certs = json.load(open("certs_min.json")); F = Field(5)
c = certs[1]; ex = certs[0]
z = int(re.search(r"delete dir (\d+)", ex["src"]).group(1))
dU = [ex["moves"][i] for i in range(len(ex["moves"])) if ex["cert"]["U"] >> i & 1] + [z]
cfg = Config(F, c["Ts"], c["moves"]); used=set(); U=[]
for d in dU:
    i = next(i for i in range(cfg.M) if cfg.moves[i]==d and i not in used); used.add(i); U.append(i)
U = sorted(U); sup = supports_in(cfg, U); R = Ring(cfg.W, U); ords = list(orders(cfg, U))
print(c["tag"], "U", U, "dirs", [cfg.moves[i] for i in U], "supports", len(sup), "orders", len(ords))
os.makedirs("jobs/tpl", exist_ok=True)
sy = [system(cfg, R, o, sup) for o in ords[:1]]
open("jobs/tpl/pl_Q1.sing","w").write(singular_text(R, sy, 0, "pl"))
open("jobs/tpl/pl_p1.sing","w").write(singular_text(R, sy, 32003, "pl"))
open("jobs/tpl/pl.ms","w").write(msolve_text(sy[0], 65521))
